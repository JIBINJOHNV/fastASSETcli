arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L || !nzchar(arguments[1L])) {
  stop("Usage: Rscript install_dependencies.R R_LIBRARY_DIRECTORY")
}

library_directory <- normalizePath(
  arguments[1L],
  mustWork = FALSE
)
dir.create(library_directory, recursive = TRUE, showWarnings = FALSE)
library_directory <- normalizePath(library_directory, mustWork = TRUE)
.libPaths(unique(c(library_directory, .libPaths())))

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages(
    "data.table",
    lib = library_directory,
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("ASSET", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
      "BiocManager",
      lib = library_directory,
      repos = "https://cloud.r-project.org"
    )
  }
  BiocManager::install(
    "ASSET",
    lib = library_directory,
    ask = FALSE,
    update = FALSE
  )
}

missing <- c("data.table", "ASSET")[!vapply(
  c("data.table", "ASSET"),
  requireNamespace,
  logical(1L),
  quietly = TRUE
)]
if (length(missing) > 0L) {
  stop("Dependency installation failed: ", paste(missing, collapse = ", "))
}

cat("Dependencies available in ", library_directory, "\n", sep = "")
