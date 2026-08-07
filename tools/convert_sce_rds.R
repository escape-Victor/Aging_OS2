#!/usr/bin/env Rscript
#
# convert_sce_rds.R — Convert a SingleCellExperiment .rds to .h5ad via zellkonverter
#
# Called by the Python subprocess bridge during Route 6 (RDS/SCE conversion).
# This script reads a serialized SingleCellExperiment object and writes an
# AnnData-compatible .h5ad file, preserving all assays, rowData, colData,
# reducedDims, and metadata.
#
# Usage:
#   Rscript convert_sce_rds.R --input  /path/to/input.rds \
#                             --output /path/to/output.h5ad \
#                             [--sample-id SAMPLE123]
#
# Exit codes:
#   0 — conversion succeeded
#   1 — conversion failed (error details written to stderr)

suppressPackageStartupMessages({
  library(methods)
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Emit a JSON status line to stdout.  The Python bridge parses this to
# confirm success and to extract provenance metadata.
emit_status <- function(status, message, details = list()) {
  payload <- list(
    status    = status,                 # "ok" | "error"
    message   = message,               # human-readable summary
    details   = details,               # arbitrary named list
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  )
  # Prefer jsonlite for robust escaping, but fall back to manual JSON so the
  # script works even in a bare R installation without jsonlite.
  if (is_pkg_installed("jsonlite")) {
    cat(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = FALSE), "\n")
  } else {
    cat(build_json_manually(payload), "\n")
  }
}

# Minimal JSON builder — handles only the flat structures emitted by this
# script (strings, numbers, booleans, named lists, character vectors).
# Not a general-purpose serializer; sufficient for our provenance payloads.
build_json_manually <- function(x) {
  if (is.null(x)) return("null")
  if (is.character(x) && length(x) == 1L) {
    return(sprintf('"%s"', gsub('"', '\\"', x, fixed = TRUE)))
  }
  if (is.numeric(x) && length(x) == 1L) {
    if (is.integer(x)) return(as.character(x))
    return(as.character(x))
  }
  if (is.logical(x) && length(x) == 1L) {
    return(if (x) "true" else "false")
  }
  if (is.list(x) && !is.null(names(x))) {
    parts <- vapply(names(x), function(nm) {
      sprintf('"%s": %s', nm, build_json_manually(x[[nm]]))
    }, character(1L), USE.NAMES = FALSE)
    return(sprintf("{%s}", paste(parts, collapse = ", ")))
  }
  if (is.character(x) && length(x) > 1L) {
    parts <- vapply(x, build_json_manually, character(1L), USE.NAMES = FALSE)
    return(sprintf("[%s]", paste(parts, collapse = ", ")))
  }
  # Fallback: string representation
  sprintf('"%s"', gsub('"', '\\"', as.character(x), fixed = TRUE))
}

# Emit an error, print to stderr, and quit with code 1.
die <- function(msg, ..., call. = FALSE) {
  full_msg <- sprintf(msg, ...)
  write(full_msg, stderr())

  # Try to emit structured JSON to stdout so the bridge can parse it even on
  # failure (some callers consume stdout exclusively).
  tryCatch(
    emit_status("error", full_msg),
    error = function(e) NULL
  )

  quit(save = "default", status = 1L, runLast = FALSE)
}

# Return TRUE when the package named `pkg` is installed.
is_pkg_installed <- function(pkg) {
  nzchar(system.file(package = pkg))
}

# ---------------------------------------------------------------------------
# Argument parsing (positional-style flags via commandArgs)
# ---------------------------------------------------------------------------

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) {
    cat("Usage: Rscript convert_sce_rds.R --input PATH --output PATH [--sample-id STR]\n",
        file = stderr())
    quit(save = "default", status = 1L, runLast = FALSE)
  }

  # Simple state-machine parser for --key value pairs.
  parsed <- list(
    input     = NULL,
    output    = NULL,
    sample_id = NULL
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    switch(key,
      `--input` = {
        i <- i + 1L
        if (i > length(args)) die("Missing value for --input")
        parsed$input <- args[i]
      },
      `--output` = {
        i <- i + 1L
        if (i > length(args)) die("Missing value for --output")
        parsed$output <- args[i]
      },
      `--sample-id` = {
        i <- i + 1L
        if (i > length(args)) die("Missing value for --sample-id")
        parsed$sample_id <- args[i]
      },
      {
        # Unknown flag — warn but continue.  GEO supplementary filenames
        # sometimes leak extra tokens; being lenient here avoids brittle
        # failures while still catching real misconfigurations downstream.
        write(sprintf("WARNING: ignoring unrecognized argument '%s'", key), stderr())
      }
    )
    i <- i + 1L
  }

  # Validate required arguments.
  if (is.null(parsed$input)) {
    die("Required argument --input is missing")
  }
  if (is.null(parsed$output)) {
    die("Required argument --output is missing")
  }

  parsed
}

# ---------------------------------------------------------------------------
# Input validation — filesystem checks
# ---------------------------------------------------------------------------

validate_paths <- function(input_path, output_path) {
  if (!file.exists(input_path)) {
    die("Input file does not exist: %s", input_path)
  }

  input_norm <- normalizePath(input_path, mustWork = TRUE)

  # Ensure the output directory exists (create if missing).
  out_dir <- dirname(output_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(out_dir)) {
      die("Cannot create output directory: %s", out_dir)
    }
    write(sprintf("Created output directory: %s", out_dir), stderr())
  }

  if (normalizePath(out_dir, mustWork = FALSE) == dirname(input_norm) &&
      basename(output_path) == basename(input_path)) {
    die("Input and output paths must differ (refusing to overwrite source)")
  }

  list(input = input_norm, output = output_path)
}

# ---------------------------------------------------------------------------
# Pre-flight: check that zellkonverter (and its heavy dependency basilisk)
# are installed.  If not, emit a clear actionable message and quit.
# ---------------------------------------------------------------------------

check_dependencies <- function() {
  if (!is_pkg_installed("zellkonverter")) {
    die(paste(
      "zellkonverter is not installed.",
      "",
      "Install it with:",
      "  BiocManager::install('zellkonverter')",
      "",
      "zellkonverter depends on basilisk, which provisions a conda environment",
      "with Python, anndata, and h5py. The first call after installation will",
      "trigger an environment build (2-5 minutes); subsequent calls are fast.",
      sep = "\n"
    ))
  }

  # zellkonverter will pull in basilisk and reticulate automatically when
  # loaded, but we suppress the startup chatter.
  suppressPackageStartupMessages({
    library(zellkonverter)
  })

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Core conversion logic
# ---------------------------------------------------------------------------

convert_sce_to_h5ad <- function(input_path, output_path, sample_id = NULL) {
  # --- Load the RDS -------------------------------------------------------
  write(sprintf("Loading RDS: %s", input_path), stderr())
  obj <- tryCatch(
    readRDS(input_path),
    error = function(e) {
      die("Failed to read RDS file '%s': %s", input_path, e$message)
    }
  )

  # --- Verify type --------------------------------------------------------
  if (!inherits(obj, "SingleCellExperiment")) {
    actual <- paste(class(obj), collapse = ", ")
    die(
      "Expected a SingleCellExperiment object but got: %s\n  File: %s\n  Classes: %s",
      actual, input_path, actual
    )
  }

  write(sprintf(
    "Object OK: SingleCellExperiment with %d genes x %d cells",
    nrow(obj), ncol(obj)
  ), stderr())

  # --- Inject / validate sample_id ----------------------------------------
  # Route 6 requires sample-level traceability.  If the caller supplied
  # --sample-id we stamp it into colData; if the object already has it via
  # colData we preserve it; if neither, we emit a warning but proceed.
  if (!is.null(sample_id) && nchar(sample_id) > 0L) {
    # Overwrite or create the column.
    SummarizedExperiment::colData(obj)$sample_id <- sample_id
    write(sprintf("Set sample_id in colData: %s", sample_id), stderr())
  } else if ("sample_id" %in% colnames(SummarizedExperiment::colData(obj))) {
    n_unique <- length(unique(SummarizedExperiment::colData(obj)$sample_id))
    write(sprintf(
      "Found existing sample_id in colData (%d unique value(s))", n_unique
    ), stderr())
  } else {
    write("WARNING: no --sample-id provided and colData has no 'sample_id' column",
          stderr())
  }

  # --- Convert ------------------------------------------------------------
  write(sprintf("Converting to H5AD: %s", output_path), stderr())
  t_start <- Sys.time()

  tryCatch(
    {
      zellkonverter::writeH5AD(sce = obj, file = output_path)
    },
    error = function(e) {
      die("zellkonverter::writeH5AD() failed: %s", e$message)
    }
  )

  t_elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 2)
  write(sprintf("Conversion completed in %.2f seconds", t_elapsed), stderr())

  # --- Verify output exists and has non-zero size -------------------------
  if (!file.exists(output_path)) {
    die("Output file was not created: %s", output_path)
  }

  out_size <- file.info(output_path)$size
  if (is.na(out_size) || out_size == 0L) {
    die("Output file is empty: %s", output_path)
  }

  write(sprintf("Output file size: %s bytes", format(out_size, big.mark = ",")),
        stderr())

  # --- Return provenance metadata for the success report ------------------
  list(
    input_class     = class(obj)[1L],
    n_genes         = nrow(obj),
    n_cells         = ncol(obj),
    assay_names     = as.character(SummarizedExperiment::assayNames(obj)),
    colData_columns = colnames(SummarizedExperiment::colData(obj)),
    reducedDim_names = names(SingleCellExperiment::reducedDims(obj)),
    sample_id       = if (!is.null(sample_id)) sample_id else
                       if ("sample_id" %in% colnames(SummarizedExperiment::colData(obj)))
                         unique(SummarizedExperiment::colData(obj)$sample_id) else NULL,
    output_size_bytes = out_size,
    elapsed_seconds   = t_elapsed
  )
}

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

main <- function() {
  # 1. Parse CLI arguments
  cli <- parse_args()
  paths <- validate_paths(cli$input, cli$output)

  # 2. Check dependencies (early fail with a clear message)
  check_dependencies()

  # 3. Convert
  provenance <- convert_sce_to_h5ad(
    input_path  = paths$input,
    output_path = paths$output,
    sample_id   = cli$sample_id
  )

  # 4. Report success (JSON to stdout)
  emit_status("ok", "Successfully converted SingleCellExperiment to H5AD", provenance)
}

# --- Run -------------------------------------------------------------------
main()
quit(save = "default", status = 0L, runLast = FALSE)
