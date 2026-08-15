make_run_signature <- function(params, trait_names, correlation_matrix) {
  input_info <- file.info(params$fastasset_input)
  ldsc_info <- file.info(params$ldsc_rdata)
  signature_values <- c(
    "pipeline_version=2026-08-15-reference-aligned-cli-meta-v2",
    paste0("input=", params$fastasset_input),
    paste0("input_size=", input_info$size),
    paste0("input_mtime=", as.numeric(input_info$mtime)),
    paste0("ldsc=", params$ldsc_rdata),
    paste0("ldsc_size=", ldsc_info$size),
    paste0("ldsc_mtime=", as.numeric(ldsc_info$mtime)),
    paste0("ldsc_object=", params$ldsc_object_name),
    paste0("traits=", paste(trait_names, collapse = ";")),
    paste0("scr_pthr=", format(params$scr_pthr, digits = 16)),
    paste0("max_side=", params$max_numtraits_per_side),
    paste0("min_traits=", params$min_available_traits),
    paste0("cor_thr=", format(params$cor_thr, digits = 16)),
    paste0("meth_pval=", params$meth_pval),
    paste0("include_meta=", params$include_meta),
    paste0("chunk_size=", params$chunk_size),
    paste0("cor=", paste(format(correlation_matrix, digits = 12), collapse = ";"))
  )

  signature_file <- tempfile("fastasset_signature_", fileext = ".txt")
  on.exit(unlink(signature_file, force = TRUE), add = TRUE)
  writeLines(signature_values, signature_file, useBytes = TRUE)
  unname(substr(tools::md5sum(signature_file), 1L, 12L))
}

order_chunk_files <- function(files) {
  if (length(files) == 0L) return(files)
  chunk_id <- suppressWarnings(
    as.integer(sub(".*chunk_([0-9]+).*", "\\1", basename(files)))
  )
  files[order(chunk_id, na.last = TRUE)]
}

combine_tsv_gz <- function(files, output_file) {
  files <- order_chunk_files(files[file.exists(files)])
  if (length(files) == 0L) return(NULL)
  if (!nzchar(Sys.which("gzip"))) stop("The gzip executable is required.")

  temporary_tsv <- paste0(output_file, ".building.tsv")
  temporary_gz <- paste0(temporary_tsv, ".gz")
  unlink(c(temporary_tsv, temporary_gz, output_file), force = TRUE)

  wrote_header <- FALSE
  for (file in files) {
    table <- data.table::fread(file)
    if (nrow(table) == 0L && wrote_header) next

    data.table::fwrite(
      table,
      temporary_tsv,
      sep = "\t",
      append = wrote_header,
      col.names = !wrote_header,
      quote = FALSE,
      na = "NA"
    )
    wrote_header <- TRUE
  }

  if (!wrote_header) return(NULL)

  gzip_status <- system2("gzip", c("-f", shQuote(temporary_tsv)))
  if (!identical(gzip_status, 0L)) {
    stop("gzip failed while creating ", output_file)
  }
  if (!file.rename(temporary_gz, output_file)) {
    stop("Could not move combined output into place: ", output_file)
  }

  output_file
}

capture_environment <- function(output_dir) {
  writeLines(
    capture.output(sessionInfo()),
    file.path(output_dir, "sessionInfo.txt")
  )

  external_tools <- c("gzip", "python", "plink", "bcftools")
  versions <- lapply(external_tools, function(tool) {
    executable <- Sys.which(tool)
    if (!nzchar(executable)) return(NA_character_)
    tryCatch(
      system2(executable, "--version", stdout = TRUE, stderr = TRUE)[1L],
      error = function(e) NA_character_
    )
  })
  tool_table <- data.table::data.table(
    tool = external_tools,
    version = unlist(versions, use.names = FALSE)
  )
  data.table::fwrite(
    tool_table,
    file.path(output_dir, "external_tools.tsv"),
    sep = "\t"
  )
  invisible(TRUE)
}
