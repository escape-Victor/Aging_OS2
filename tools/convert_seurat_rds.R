#!/usr/bin/env Rscript
#
# convert_seurat_rds.R — Convert Seurat or SingleCellExperiment .rds to .h5ad
#
# Called by the Python subprocess bridge in convert_rds_to_h5ad().
# Dispatched by: micromamba run -n r-converters Rscript tools/convert_seurat_rds.R
#
# Usage:
#   Rscript convert_seurat_rds.R --input <path> --output <path> \
#       [--sample-id <str>] [--assay <str>]
#
# Exit codes:
#   0 — success (JSON status line written to stdout)
#   1 — error   (message written to stderr)
# ============================================================================

suppressPackageStartupMessages({
  library(methods)   # always available; provides is(), class()
})

# ---------------------------------------------------------------------------
# 1. Parse command-line arguments
# ---------------------------------------------------------------------------

argv <- commandArgs(trailingOnly = TRUE)

# Simple argument parser — handles --key value (no external deps)
parse_args <- function(args) {
  # Normalise: split "--key=value" into "--key" "value"
  expanded <- character(0)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (grepl("^--[^=]+=", arg)) {
      parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
      expanded <- c(expanded, parts[[1L]], parts[[2L]])
    } else {
      expanded <- c(expanded, arg)
    }
    i <- i + 1L
  }

  parsed <- list()
  i <- 1L
  while (i <= length(expanded)) {
    key <- expanded[[i]]
    if (grepl("^--", key)) {
      key_name <- sub("^--", "", key)
      next_idx <- i + 1L
      if (next_idx <= length(expanded) && !grepl("^--", expanded[[next_idx]])) {
        parsed[[key_name]] <- expanded[[next_idx]]
        i <- next_idx
      } else {
        # Boolean flag (no value follows)
        parsed[[key_name]] <- TRUE
      }
    }
    i <- i + 1L
  }
  parsed
}

opts <- parse_args(argv)

# --- Validate required arguments ---
if (is.null(opts[["input"]])) {
  write("ERROR: --input is required (path to .rds file)", stderr())
  quit(status = 1L)
}
if (is.null(opts[["output"]])) {
  write("ERROR: --output is required (path for output .h5ad file)", stderr())
  quit(status = 1L)
}

input_path  <- opts[["input"]]
output_path <- opts[["output"]]
sample_id   <- opts[["sample-id"]]    # NULL if omitted
assay_name  <- opts[["assay"]]        # NULL if omitted
if (is.null(assay_name)) assay_name <- "RNA"

# Validate input file exists
if (!file.exists(input_path)) {
  write(sprintf("ERROR: input file not found: %s", input_path), stderr())
  quit(status = 1L)
}

# Ensure output directory exists
output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# ---------------------------------------------------------------------------
# 2. Load required packages — fail early with clear messages
# ---------------------------------------------------------------------------

load_pkg <- function(pkg_name, required_version = NULL) {
  if (!suppressPackageStartupMessages(requireNamespace(pkg_name, quietly = TRUE))) {
    write(sprintf(
      "ERROR: Package '%s' is not installed. Install it with:\n  install.packages(\"%s\")",
      pkg_name, pkg_name
    ), stderr())
    quit(status = 1L)
  }
  suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
  if (!is.null(required_version)) {
    ver <- as.character(packageVersion(pkg_name))
    if (utils::compareVersion(ver, required_version) < 0L) {
      write(sprintf(
        "ERROR: Package '%s' version %s is too old (need >= %s). Run:\n  install.packages(\"%s\")",
        pkg_name, ver, required_version, pkg_name
      ), stderr())
      quit(status = 1L)
    }
  }
  invisible(NULL)
}

# --- Load Seurat ecosystem packages ---
# SeuratDisk > 0.0.0.9019 fixes critical bugs in h5ad conversion
load_pkg("SeuratObject")
load_pkg("Seurat")
load_pkg("SeuratDisk", required_version = "0.0.0.9019")

# --- Lazy-load SCE support packages (only loaded if needed) ---
ensure_zellkonverter <- function() {
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    write(
      "ERROR: Package 'zellkonverter' is not installed. Install it with:\n  BiocManager::install(\"zellkonverter\")",
      stderr()
    )
    quit(status = 1L)
  }
  suppressPackageStartupMessages(library("zellkonverter", character.only = TRUE))
}

# ---------------------------------------------------------------------------
# 3. Check SeuratDisk version compatibility
# ---------------------------------------------------------------------------

check_seuratdisk_compat <- function() {
  sd_ver <- as.character(packageVersion("SeuratDisk"))
  seurat_ver <- as.character(packageVersion("Seurat"))

  # SeuratDisk >= 0.0.0.9021 requires Seurat >= 5.1.0
  if (utils::compareVersion(sd_ver, "0.0.0.9021") >= 0L &&
      utils::compareVersion(seurat_ver, "5.1.0") < 0L) {
    write(sprintf(
      "ERROR: SeuratDisk %s requires Seurat >= 5.1.0, but Seurat %s is installed.\n  ",
      sd_ver, seurat_ver
    ), stderr())
    quit(status = 1L)
  }

  message(sprintf(
    "[convert_seurat_rds] SeuratDisk %s + Seurat %s — version check OK",
    sd_ver, seurat_ver
  ))
}

check_seuratdisk_compat()

# ---------------------------------------------------------------------------
# 4. Load the .rds file
# ---------------------------------------------------------------------------

message(sprintf("[convert_seurat_rds] Loading: %s", input_path))

obj <- tryCatch(
  readRDS(input_path),
  error = function(e) {
    write(sprintf("ERROR: Failed to readRDS('%s'): %s", input_path, e$message), stderr())
    quit(status = 1L)
  }
)

obj_class <- class(obj)

# ---------------------------------------------------------------------------
# 5. Detect object type and dispatch conversion
# ---------------------------------------------------------------------------

# --- Helper: emit JSON status line on success ---
emit_success <- function(out_path, obj_type, n_cells, n_genes) {
  status <- list(
    success       = TRUE,
    output_path   = out_path,
    object_type   = obj_type,
    n_cells       = n_cells,
    n_genes       = n_genes,
    provenance = list(
      tool     = "convert_seurat_rds.R",
      input    = input_path,
      assay    = assay_name,
      sample_id = if (is.null(sample_id)) NULL else sample_id
    )
  )
  # Build JSON manually — no external dependency
  json_parts <- character(0)
  json_parts <- c(json_parts, sprintf('  "success": %s', if (status$success) "true" else "false"))
  json_parts <- c(json_parts, sprintf('  "output_path": "%s"', status$output_path))
  json_parts <- c(json_parts, sprintf('  "object_type": "%s"', status$object_type))
  json_parts <- c(json_parts, sprintf('  "n_cells": %d', status$n_cells))
  json_parts <- c(json_parts, sprintf('  "n_genes": %d', status$n_genes))
  json_parts <- c(json_parts, sprintf('  "provenance": {'))
  json_parts <- c(json_parts, sprintf('    "tool": "%s"', status$provenance$tool))
  json_parts <- c(json_parts, sprintf('    "input": "%s"', status$provenance$input))
  json_parts <- c(json_parts, sprintf('    "assay": "%s"', status$provenance$assay))
  if (is.null(status$provenance$sample_id)) {
    json_parts <- c(json_parts, '    "sample_id": null')
  } else {
    json_parts <- c(json_parts, sprintf('    "sample_id": "%s"', status$provenance$sample_id))
  }
  json_parts <- c(json_parts, "  }")
  cat("{\n", paste(json_parts, collapse = ",\n"), "\n}\n", sep = "")
}

# ---------------------------------------------------------------------------
# 5a. Seurat object → h5ad via SeuratDisk
# ---------------------------------------------------------------------------
convert_seurat <- function(seurat_obj) {
  message(sprintf("[convert_seurat_rds] Detected Seurat object (class: %s)", paste(class(seurat_obj), collapse = ", ")))

  # Create a temporary .h5seurat file in the output directory
  # (same filesystem → efficient rename, no cross-volume issues)
  tmp_h5seurat <- file.path(output_dir, paste0(".tmp_", basename(tempfile(pattern = "seurat_")), ".h5seurat"))

  message(sprintf("[convert_seurat_rds] Saving intermediate h5Seurat: %s", tmp_h5seurat))

  # Step 1: Save as .h5seurat
  tryCatch({
    SaveH5Seurat(seurat_obj, filename = tmp_h5seurat, overwrite = TRUE)
  }, error = function(e) {
    if (file.exists(tmp_h5seurat)) file.remove(tmp_h5seurat)
    write(sprintf("ERROR: SaveH5Seurat() failed: %s", e$message), stderr())
    quit(status = 1L)
  })

  # Step 2: Convert .h5seurat → .h5ad
  # Convert() writes to the same directory as the input by default,
  # appending ".h5ad" extension. We rename to the requested output path.
  message(sprintf("[convert_seurat_rds] Converting h5Seurat → h5ad"))

  tryCatch({
    Convert(tmp_h5seurat, dest = "h5ad", overwrite = TRUE)
  }, error = function(e) {
    # Clean up temp files
    if (file.exists(tmp_h5seurat)) file.remove(tmp_h5seurat)
    default_h5ad <- paste0(tmp_h5seurat, ".h5ad")
    if (file.exists(default_h5ad)) file.remove(default_h5ad)
    write(sprintf("ERROR: Convert() h5seurat→h5ad failed: %s", e$message), stderr())
    quit(status = 1L)
  })

  # Step 3: Rename to target output path
  converted_h5ad <- paste0(tmp_h5seurat, ".h5ad")
  if (!file.exists(converted_h5ad)) {
    if (file.exists(tmp_h5seurat)) file.remove(tmp_h5seurat)
    write(sprintf("ERROR: Convert() reported success but output file not found: %s", converted_h5ad), stderr())
    quit(status = 1L)
  }

  success <- file.rename(converted_h5ad, output_path)
  # Clean up temp .h5seurat (Convert leaves it behind)
  if (file.exists(tmp_h5seurat)) file.remove(tmp_h5seurat)

  if (!success) {
    # Fallback: copy + remove
    if (!file.copy(converted_h5ad, output_path, overwrite = TRUE)) {
      write(sprintf("ERROR: Failed to move '%s' → '%s'", converted_h5ad, output_path), stderr())
      quit(status = 1L)
    }
    file.remove(converted_h5ad)
  }

  # Return dimensions for reporting
  list(
    n_cells = ncol(seurat_obj),
    n_genes = nrow(seurat_obj)
  )
}

# ---------------------------------------------------------------------------
# 5b. SingleCellExperiment → h5ad via zellkonverter
# ---------------------------------------------------------------------------
convert_sce <- function(sce_obj) {
  ensure_zellkonverter()

  message(sprintf("[convert_seurat_rds] Detected SingleCellExperiment (class: %s)", paste(class(sce_obj), collapse = ", ")))

  tryCatch({
    zellkonverter::writeH5AD(
      sce_obj,
      file    = output_path,
      assay_name = assay_name
    )
  }, error = function(e) {
    write(sprintf("ERROR: zellkonverter::writeH5AD() failed: %s", e$message), stderr())
    quit(status = 1L)
  })

  list(
    n_cells = ncol(sce_obj),
    n_genes = nrow(sce_obj)
  )
}

# ---------------------------------------------------------------------------
# 6. Dispatch based on detected class
# ---------------------------------------------------------------------------

is_seurat <- "Seurat" %in% obj_class
is_sce    <- "SingleCellExperiment" %in% obj_class

if (!is_seurat && !is_sce) {
  # Attempt to infer: some packages add their own subclass
  # e.g. "Seurat" may appear in inherits() even if not direct class member
  if (inherits(obj, "Seurat")) {
    is_seurat <- TRUE
  } else if (inherits(obj, "SingleCellExperiment")) {
    is_sce <- TRUE
  } else {
    write(sprintf(
      "ERROR: Unknown R object type. Expected 'Seurat' or 'SingleCellExperiment', got: %s",
      paste(obj_class, collapse = ", ")
    ), stderr())
    quit(status = 1L)
  }
}

if (is_seurat) {
  dims <- convert_seurat(obj)
} else {
  dims <- convert_sce(obj)
}

# ---------------------------------------------------------------------------
# 7. Verify output and emit success
# ---------------------------------------------------------------------------

if (!file.exists(output_path)) {
  write(sprintf("ERROR: Conversion appeared to succeed but output file is missing: %s", output_path), stderr())
  quit(status = 1L)
}

output_size <- file.info(output_path)$size
message(sprintf(
  "[convert_seurat_rds] Success: %s (%d cells x %d genes, %.1f MB)",
  output_path, dims$n_cells, dims$n_genes, output_size / 1e6
))

emit_success(output_path, if (is_seurat) "seurat" else "sce", dims$n_cells, dims$n_genes)
quit(status = 0L)
