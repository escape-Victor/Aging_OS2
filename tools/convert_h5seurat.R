#!/usr/bin/env Rscript
#
# convert_h5seurat.R — Convert .h5seurat to .h5ad via SeuratDisk
#
# Called by the Python subprocess bridge (r-converters environment).
# Reads a Seurat v5 .h5seurat file and writes an AnnData .h5ad file.
#
# Usage:
#   Rscript convert_h5seurat.R --input <path> --output <path> [--sample-id <id>]
#
# Exit codes:
#   0 — success (JSON status written to stdout)
#   1 — error (message written to stderr)

suppressPackageStartupMessages({
  library(SeuratDisk)
})

# ---------------------------------------------------------------------------
# 1. Parse command-line arguments
# ---------------------------------------------------------------------------

parse_args <- function(args) {
  # args: character vector from commandArgs(trailingOnly = TRUE)

  if (length(args) == 0L) {
    stop("No arguments provided. Use --input PATH --output PATH [--sample-id STR]")
  }

  input_path  <- NULL
  output_path <- NULL
  sample_id   <- NULL

  i <- 1L
  while (i <= length(args)) {
    switch(
      args[i],
      "--input" = {
        if (i + 1L > length(args)) stop("--input requires a path argument")
        input_path <- args[i + 1L]
        i <- i + 1L
      },
      "--output" = {
        if (i + 1L > length(args)) stop("--output requires a path argument")
        output_path <- args[i + 1L]
        i <- i + 1L
      },
      "--sample-id" = {
        if (i + 1L > length(args)) stop("--sample-id requires a string argument")
        sample_id <- args[i + 1L]
        i <- i + 1L
      },
      {
        stop(sprintf("Unknown argument: %s", args[i]))
      }
    )
    i <- i + 1L
  }

  list(input = input_path, output = output_path, sample_id = sample_id)
}

# ---------------------------------------------------------------------------
# 2. Input validation
# ---------------------------------------------------------------------------

validate_inputs <- function(parsed) {
  # Check required arguments
  if (is.null(parsed$input) || nchar(parsed$input) == 0L) {
    stop("Missing required argument: --input PATH")
  }
  if (is.null(parsed$output) || nchar(parsed$output) == 0L) {
    stop("Missing required argument: --output PATH")
  }

  # Check that input file exists
  if (!file.exists(parsed$input)) {
    stop(sprintf("Input file does not exist: %s", parsed$input))
  }

  # Check that input has .h5seurat extension
  if (!grepl("\\.h5seurat$", parsed$input, ignore.case = TRUE)) {
    stop(sprintf(
      "Input file must have .h5seurat extension, got: %s",
      basename(parsed$input)
    ))
  }

  # Check that output has .h5ad extension
  if (!grepl("\\.h5ad$", parsed$output, ignore.case = TRUE)) {
    stop(sprintf(
      "Output file must have .h5ad extension, got: %s",
      basename(parsed$output)
    ))
  }

  # Warn if input and output are the same path
  if (normalizePath(parsed$input, mustWork = FALSE) ==
      normalizePath(parsed$output, mustWork = FALSE)) {
    stop("Input and output paths must be different")
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# 3. Convert via SeuratDisk
# ---------------------------------------------------------------------------

convert_h5seurat_to_h5ad <- function(input_path, output_path, sample_id = NULL) {
  message(sprintf("Converting %s -> %s", input_path, output_path))

  # SeuratDisk::Convert handles both directions. When dest is a path ending
  # in .h5ad, it writes directly to that path. When dest is "h5ad" it
  # auto-generates the filename by replacing the extension.
  #
  # We pass the full output path as dest to control the output location.
  SeuratDisk::Convert(
    source    = input_path,
    dest      = output_path,
    overwrite = TRUE
  )

  # Verify the output was created
  if (!file.exists(output_path)) {
    stop(sprintf("Conversion did not produce output file: %s", output_path))
  }

  # Collect output metadata
  output_size <- file.info(output_path)$size
  result <- list(
    status       = "success",
    input        = normalizePath(input_path, mustWork = TRUE),
    output       = normalizePath(output_path, mustWork = TRUE),
    output_bytes = as.integer(output_size)
  )

  if (!is.null(sample_id) && nchar(sample_id) > 0L) {
    result$sample_id <- sample_id
  }

  # Write JSON status line to stdout for the Python bridge to consume
  cat(jsonlite::toJSON(result, auto_unbox = TRUE), "\n", sep = "")

  invisible(result)
}

# ---------------------------------------------------------------------------
# 4. Main entry point
# ---------------------------------------------------------------------------

main <- function() {
  # Capture all errors and emit structured output
  tryCatch(
    {
      args   <- parse_args(commandArgs(trailingOnly = TRUE))
      validate_inputs(args)
      convert_h5seurat_to_h5ad(args$input, args$output, args$sample_id)
      quit(status = 0L, save = "no")
    },
    error = function(e) {
      # Write error to stderr so the Python bridge can capture it
      err_json <- jsonlite::toJSON(
        list(status = "error", message = conditionMessage(e)),
        auto_unbox = TRUE
      )
      writeLines(err_json, con = stderr())
      quit(status = 1L, save = "no")
    }
  )
}

main()
