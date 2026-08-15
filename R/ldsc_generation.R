normalize_manifest_header <- function(values) {
  normalized <- toupper(gsub("[^A-Za-z0-9]+", "_", trimws(values)))
  gsub("^_+|_+$", "", normalized)
}

find_manifest_column <- function(headers, aliases, label, required = FALSE) {
  hits <- which(headers %in% aliases)
  if (length(hits) > 1L) {
    stop(
      "The LDSC manifest has multiple columns that map to ", label, ": ",
      paste(headers[hits], collapse = ", "), "."
    )
  }
  if (length(hits) == 0L) {
    if (required) {
      stop(
        "The LDSC manifest is missing ", label, ". Accepted header(s): ",
        paste(aliases, collapse = ", "), "."
      )
    }
    return(NA_integer_)
  }
  hits
}

coerce_manifest_number <- function(values, label) {
  text <- trimws(as.character(values))
  missing <- is.na(values) | !nzchar(text) |
    toupper(text) %in% c("NA", "N/A", "NULL", ".")
  parsed <- suppressWarnings(as.numeric(text))
  invalid <- !missing & !is.finite(parsed)
  if (any(invalid)) {
    stop(
      "The LDSC manifest has invalid ", label, " value(s) on row(s): ",
      paste(which(invalid) + 1L, collapse = ", "), "."
    )
  }
  parsed[missing] <- NA_real_
  parsed
}

resolve_manifest_paths <- function(paths, manifest_path) {
  paths <- path.expand(trimws(as.character(paths)))
  if (anyNA(paths) || any(!nzchar(paths))) {
    stop("Every LDSC manifest FILE value must be non-empty.")
  }

  is_absolute <- grepl("^/", paths) | grepl("^[A-Za-z]:[/\\\\]", paths)
  paths[!is_absolute] <- file.path(dirname(manifest_path), paths[!is_absolute])
  missing <- !file.exists(paths)
  if (any(missing)) {
    stop(
      "Summary-statistic file(s) listed in the LDSC manifest do not exist:\n  ",
      paste(paths[missing], collapse = "\n  ")
    )
  }
  unreadable <- file.access(paths, 4L) != 0L
  if (any(unreadable)) {
    stop(
      "Summary-statistic file(s) listed in the LDSC manifest are not readable:\n  ",
      paste(paths[unreadable], collapse = "\n  ")
    )
  }
  vapply(paths, normalizePath, character(1L), mustWork = TRUE)
}

read_ldsc_manifest <- function(path, analysis_traits) {
  analysis_traits <- as.character(analysis_traits)
  manifest <- data.table::fread(
    path,
    na.strings = c("", "NA", "N/A", "NULL", "."),
    check.names = FALSE
  )
  if (nrow(manifest) == 0L) stop("The LDSC manifest has no trait rows.")

  headers <- normalize_manifest_header(names(manifest))
  if (anyDuplicated(headers)) {
    stop("The LDSC manifest has duplicate normalized column names.")
  }
  trait_column <- find_manifest_column(
    headers, c("TRAIT", "NAME", "TRAIT_NAME"), "TRAIT", required = TRUE
  )
  file_column <- find_manifest_column(
    headers, c("FILE", "PATH", "SUMSTATS", "SUMSTATS_FILE"),
    "FILE", required = TRUE
  )
  n_column <- find_manifest_column(
    headers, c("N", "SAMPLE_SIZE", "TOTAL_N"), "N"
  )
  sample_prev_column <- find_manifest_column(
    headers, c("SAMPLE_PREV", "SPREV", "SAMPLE_PREVALENCE"), "SAMPLE_PREV"
  )
  population_prev_column <- find_manifest_column(
    headers,
    c("POPULATION_PREV", "PPREV", "POPULATION_PREVALENCE"),
    "POPULATION_PREV"
  )

  traits <- trimws(as.character(manifest[[trait_column]]))
  if (anyNA(traits) || any(!nzchar(traits))) {
    stop("Every LDSC manifest TRAIT value must be non-empty.")
  }
  if (anyDuplicated(traits)) {
    stop(
      "Duplicate LDSC manifest trait(s): ",
      paste(unique(traits[duplicated(traits)]), collapse = ", "), "."
    )
  }

  missing_traits <- setdiff(analysis_traits, traits)
  extra_traits <- setdiff(traits, analysis_traits)
  trait_errors <- c()
  if (length(missing_traits) > 0L) {
    trait_errors <- c(
      trait_errors,
      paste0("missing from manifest: ", paste(missing_traits, collapse = ", "))
    )
  }
  if (length(extra_traits) > 0L) {
    trait_errors <- c(
      trait_errors,
      paste0("not present in FastASSET input: ", paste(extra_traits, collapse = ", "))
    )
  }
  if (length(trait_errors) > 0L) {
    stop(
      "Automatic LDSC requires exactly one manifest row for every <trait>.Beta ",
      "column:\n  ", paste(trait_errors, collapse = "\n  ")
    )
  }

  input_files <- resolve_manifest_paths(manifest[[file_column]], path)
  if (anyDuplicated(input_files)) {
    stop("Each LDSC manifest trait must use a distinct summary-statistic file.")
  }

  sample_size <- if (is.na(n_column)) {
    rep(NA_real_, nrow(manifest))
  } else {
    coerce_manifest_number(manifest[[n_column]], "N")
  }
  if (any(!is.na(sample_size) & sample_size <= 0)) {
    stop("Every non-missing LDSC manifest N value must be greater than zero.")
  }

  sample_prev <- if (is.na(sample_prev_column)) {
    rep(NA_real_, nrow(manifest))
  } else {
    coerce_manifest_number(manifest[[sample_prev_column]], "SAMPLE_PREV")
  }
  population_prev <- if (is.na(population_prev_column)) {
    rep(NA_real_, nrow(manifest))
  } else {
    coerce_manifest_number(
      manifest[[population_prev_column]], "POPULATION_PREV"
    )
  }
  unmatched_prevalence <- xor(is.na(sample_prev), is.na(population_prev))
  if (any(unmatched_prevalence)) {
    stop(
      "SAMPLE_PREV and POPULATION_PREV must both be NA for quantitative ",
      "traits or both be provided for binary traits. Invalid manifest row(s): ",
      paste(which(unmatched_prevalence) + 1L, collapse = ", "), "."
    )
  }
  invalid_prevalence <-
    (!is.na(sample_prev) & (sample_prev <= 0 | sample_prev >= 1)) |
    (!is.na(population_prev) &
       (population_prev <= 0 | population_prev >= 1))
  if (any(invalid_prevalence)) {
    stop("Every supplied prevalence must be strictly between zero and one.")
  }

  resolved <- data.table::data.table(
    trait = traits,
    file = input_files,
    N = sample_size,
    sample_prev = sample_prev,
    population_prev = population_prev
  )
  resolved[["order"]] <- match(resolved[["trait"]], analysis_traits)
  if (anyNA(resolved[["order"]])) {
    stop("Internal error while reordering the LDSC manifest.")
  }
  data.table::setorder(resolved, order)
  data.table::setcolorder(
    resolved,
    c("order", "trait", "file", "N", "sample_prev", "population_prev")
  )
  resolved
}

ld_reference_paths <- function(directory, chromosomes, include_m = FALSE) {
  score_files <- file.path(
    directory,
    paste0(seq_len(chromosomes), ".l2.ldscore.gz")
  )
  if (!include_m) return(score_files)
  c(
    score_files,
    file.path(directory, paste0(seq_len(chromosomes), ".l2.M_5_50"))
  )
}

validate_ld_reference <- function(params) {
  ld_files <- ld_reference_paths(
    params$ld_ref, params$ldsc_chr, include_m = TRUE
  )
  separate_weights <- !identical(params$ld_ref, params$wld_ref)
  weight_files <- if (separate_weights) {
    ld_reference_paths(params$wld_ref, params$ldsc_chr, include_m = FALSE)
  } else {
    character()
  }
  required <- c(ld_files, weight_files)
  missing <- !file.exists(required)
  if (any(missing)) {
    stop(
      "Required GenomicSEM LD-score reference file(s) are missing:\n  ",
      paste(required[missing], collapse = "\n  ")
    )
  }
  unreadable <- file.access(required, 4L) != 0L
  if (any(unreadable)) {
    stop(
      "GenomicSEM LD-score reference file(s) are not readable:\n  ",
      paste(required[unreadable], collapse = "\n  ")
    )
  }
  required
}

file_metadata_lines <- function(label, paths) {
  information <- file.info(paths)
  paste(
    label,
    normalizePath(paths, mustWork = TRUE),
    information$size,
    format(as.numeric(information$mtime), scientific = FALSE),
    sep = "="
  )
}

make_ldsc_generation_signature <- function(params, manifest, reference_files) {
  mapping <- if (length(params$munge_column_map) == 0L) {
    "AUTO"
  } else {
    paste(
      names(params$munge_column_map),
      unlist(params$munge_column_map, use.names = FALSE),
      sep = "=",
      collapse = ","
    )
  }
  values <- c(
    "ldsc_generation_version=2026-08-15-v1",
    paste0("traits=", paste(manifest$trait, collapse = ";")),
    paste0("N=", paste(manifest$N, collapse = ";")),
    paste0("sample_prev=", paste(manifest$sample_prev, collapse = ";")),
    paste0(
      "population_prev=", paste(manifest$population_prev, collapse = ";")
    ),
    paste0("hm3=", params$hm3),
    paste0("ld_ref=", params$ld_ref),
    paste0("wld_ref=", params$wld_ref),
    paste0("info_filter=", format(params$munge_info_filter, digits = 16)),
    paste0("maf_filter=", format(params$munge_maf_filter, digits = 16)),
    paste0("column_map=", mapping),
    paste0("chr=", params$ldsc_chr),
    paste0("blocks=", params$ldsc_blocks),
    paste0("chisq_max=", params$ldsc_chisq_max),
    file_metadata_lines("manifest", params$sumstats_manifest),
    file_metadata_lines("hm3", params$hm3),
    file_metadata_lines("summary", manifest$file),
    file_metadata_lines("reference", reference_files)
  )
  make_text_signature(values)
}

add_ldsc_dimnames <- function(object, traits) {
  if (!is.list(object) || is.null(object$I)) {
    stop("GenomicSEM::ldsc() did not return a list containing matrix I.")
  }
  intercept <- as.matrix(object$I)
  if (nrow(intercept) != length(traits) ||
      ncol(intercept) != length(traits)) {
    stop("GenomicSEM::ldsc() returned I with unexpected dimensions.")
  }
  dimnames(intercept) <- list(traits, traits)
  object$I <- intercept

  for (matrix_name in c("S", "S_Stand")) {
    value <- object[[matrix_name]]
    if (!is.null(value) && is.matrix(value) &&
        nrow(value) == length(traits) && ncol(value) == length(traits)) {
      dimnames(value) <- list(traits, traits)
      object[[matrix_name]] <- value
    }
  }
  object
}

save_generated_ldsc <- function(object, object_name, path) {
  environment <- new.env(parent = emptyenv())
  assign(object_name, object, envir = environment)
  temporary <- file.path(
    dirname(path),
    paste0(".", basename(path), ".building-", Sys.getpid())
  )
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  save(
    list = object_name,
    file = temporary,
    envir = environment,
    compress = "gzip"
  )
  if (file.exists(path)) {
    backup <- paste0(
      path, ".incomplete-", format(Sys.time(), "%Y%m%dT%H%M%S")
    )
    if (!file.rename(path, backup)) {
      stop("Could not preserve the pre-existing incomplete LDSC object: ", path)
    }
    warning("Preserved a pre-existing incomplete LDSC object as: ", backup)
  }
  if (!file.rename(temporary, path)) {
    stop("Could not move the generated LDSC object into place: ", path)
  }
  invisible(path)
}

generate_ldsc_object <- function(params, analysis_traits) {
  if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
    stop(
      "Automatic-LDSC mode requires the GenomicSEM package. Re-run install.sh ",
      "without FASTASSET_SKIP_DEPENDENCIES=1 or install GenomicSEM first."
    )
  }

  manifest <- read_ldsc_manifest(params$sumstats_manifest, analysis_traits)
  reference_files <- validate_ld_reference(params)
  signature <- make_ldsc_generation_signature(
    params, manifest, reference_files
  )
  generation_dir <- file.path(
    params$output_dir,
    paste0(params$run_name, "_generated_ldsc_", signature)
  )
  munge_dir <- file.path(generation_dir, "munged")
  dir.create(munge_dir, recursive = TRUE, showWarnings = FALSE)

  manifest[, munged_prefix := file.path(
    munge_dir, sprintf("trait_%04d", order)
  )]
  manifest[, munged_file := paste0(munged_prefix, ".sumstats.gz")]
  resolved_manifest_file <- file.path(
    generation_dir, "ldsc_manifest_resolved.tsv"
  )
  data.table::fwrite(
    manifest,
    resolved_manifest_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  rdata <- file.path(
    generation_dir,
    paste0(params$run_name, "_LDSCoutput.RData")
  )
  completion_file <- file.path(generation_dir, "generated_ldsc.complete")
  complete <- file.exists(rdata) && file.exists(completion_file) &&
    identical(trimws(readLines(completion_file, warn = FALSE)[1L]), signature)
  if (complete) {
    message("Validated generated LDSC output exists; reusing: ", rdata)
    return(list(
      path = normalizePath(rdata, mustWork = TRUE),
      generation_dir = generation_dir,
      signature = signature,
      reused = TRUE
    ))
  }

  needs_munge <- !file.exists(manifest$munged_file)
  if (any(needs_munge)) {
    rows <- which(needs_munge)
    message(
      "Munging ", length(rows), " of ", nrow(manifest),
      " summary-statistic files with ",
      min(params$munge_cores, length(rows)), " worker(s)."
    )
    GenomicSEM::munge(
      files = manifest$file[rows],
      hm3 = params$hm3,
      trait.names = manifest$munged_prefix[rows],
      N = manifest$N[rows],
      info.filter = params$munge_info_filter,
      maf.filter = params$munge_maf_filter,
      log.name = file.path(generation_dir, paste0(params$run_name, "_munge")),
      column.names = params$munge_column_map,
      parallel = length(rows) > 1L && params$munge_cores > 1L,
      cores = min(params$munge_cores, length(rows)),
      overwrite = FALSE
    )
  } else {
    message("All signature-matched munged files exist; skipping munging.")
  }
  missing_munged <- !file.exists(manifest$munged_file)
  if (any(missing_munged)) {
    stop(
      "GenomicSEM::munge() did not create expected file(s):\n  ",
      paste(manifest$munged_file[missing_munged], collapse = "\n  ")
    )
  }

  number_traits <- nrow(manifest)
  effective_blocks <- if (number_traits > 18L) {
    ((number_traits + 1L) * (number_traits + 2L)) / 2 + 1L
  } else {
    params$ldsc_blocks
  }
  if (number_traits > 18L) {
    warning(
      "Current GenomicSEM overrides --ldsc-blocks for more than 18 traits. ",
      "For ", number_traits, " traits it will use ", effective_blocks,
      " jackknife blocks. Munging is parallel, but GenomicSEM::ldsc() itself ",
      "does not expose a worker-count argument."
    )
  }

  message(
    "Running GenomicSEM LDSC for ", number_traits,
    " traits in exact FastASSET input order."
  )
  LDSCoutput <- GenomicSEM::ldsc(
    traits = manifest$munged_file,
    sample.prev = manifest$sample_prev,
    population.prev = manifest$population_prev,
    ld = params$ld_ref,
    wld = params$wld_ref,
    trait.names = manifest$trait,
    sep_weights = !identical(params$ld_ref, params$wld_ref),
    chr = params$ldsc_chr,
    n.blocks = params$ldsc_blocks,
    ldsc.log = file.path(generation_dir, params$run_name),
    stand = FALSE,
    chisq.max = params$ldsc_chisq_max
  )
  LDSCoutput <- add_ldsc_dimnames(LDSCoutput, manifest$trait)
  save_generated_ldsc(LDSCoutput, params$ldsc_object_name, rdata)

  provenance <- data.table::data.table(
    parameter = c(
      "generation_signature", "number_traits", "requested_ldsc_blocks",
      "effective_ldsc_blocks", "ldsc_chr", "ld_ref", "wld_ref", "hm3",
      "manifest", "info_filter", "maf_filter", "munge_cores",
      "ldsc_chisq_max", "ldsc_object_name", "generated_rdata"
    ),
    value = as.character(c(
      signature, number_traits, params$ldsc_blocks, effective_blocks,
      params$ldsc_chr, params$ld_ref, params$wld_ref, params$hm3,
      params$sumstats_manifest, params$munge_info_filter,
      params$munge_maf_filter, params$munge_cores,
      if (is.na(params$ldsc_chisq_max)) "AUTO" else params$ldsc_chisq_max,
      params$ldsc_object_name, rdata
    ))
  )
  data.table::fwrite(
    provenance,
    file.path(generation_dir, "ldsc_generation_provenance.tsv"),
    sep = "\t",
    quote = FALSE
  )
  writeLines(signature, completion_file, useBytes = TRUE)
  message("Generated GenomicSEM LDSC object saved separately: ", rdata)

  list(
    path = normalizePath(rdata, mustWork = TRUE),
    generation_dir = generation_dir,
    signature = signature,
    reused = FALSE
  )
}

resolve_ldsc_source <- function(params, analysis_traits) {
  if (params$ldsc_mode == "existing") {
    return(list(
      path = params$ldsc_rdata,
      generation_dir = NA_character_,
      signature = NA_character_,
      reused = NA
    ))
  }
  generate_ldsc_object(params, analysis_traits)
}
