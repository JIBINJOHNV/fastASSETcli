#' Run the fastASSET command-line interface
#'
#' This is the programmatic entry point used by the installed `fastasset`
#' shell command. All scientific and compute settings are supplied through
#' command-line options; no configuration file is read.
#'
#' @param arguments Character vector of command-line arguments.
#' @return Invisibly returns the output-path list from the completed pipeline,
#'   or `NULL` after displaying help or version information.
#' @export
fastasset_cli_main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  if (any(arguments %in% c("-h", "--help"))) {
    cat(fastasset_cli_help())
    return(invisible(NULL))
  }
  if (any(arguments == "--version")) {
    cat("fastASSETcli ", as.character(utils::packageVersion("fastASSETcli")), "\n",
        sep = "")
    return(invisible(NULL))
  }

  params <- parse_fastasset_cli(arguments)

  # One internal math-library thread per fork avoids oversubscription when the
  # requested 90 chunk workers are active.
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )

  outputs <- run_asset_pipeline(params)
  cat("\nCompleted outputs:\n")
  cat("  input:   ", outputs$analysis_input, "\n", sep = "")
  cat("  LDSC:    ", outputs$ldsc_rdata, "\n", sep = "")
  cat("  results: ", outputs$results, "\n", sep = "")
  cat("  meta:    ", outputs$meta, "\n", sep = "")
  cat("  QC:      ", outputs$qc, "\n", sep = "")
  cat("  summary: ", outputs$summary, "\n", sep = "")
  invisible(outputs)
}
