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

is_vcf_summary_path <- function(path) {
  grepl("\\.(vcf|vcf\\.gz|vcf\\.bgz|bcf)$", tolower(path))
}

coerce_manifest_text <- function(values) {
  text <- trimws(as.character(values))
  missing <- is.na(values) | !nzchar(text) |
    toupper(text) %in% c("NA", "N/A", "NULL", ".")
  text[missing] <- NA_character_
  text
}

read_ldsc_manifest <- function(path, analysis_traits = NULL) {
  if (!is.null(analysis_traits)) analysis_traits <- as.character(analysis_traits)
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
  sample_column <- find_manifest_column(
    headers, c("SAMPLE", "VCF_SAMPLE", "SAMPLE_ID"), "SAMPLE"
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

  if (!is.null(analysis_traits)) {
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
        paste0("not present in analysis input: ", paste(extra_traits, collapse = ", "))
      )
    }
    if (length(trait_errors) > 0L) {
      stop(
        "The manifest and prepared analysis input must contain exactly the ",
        "same traits:\n  ", paste(trait_errors, collapse = "\n  ")
      )
    }
  }

  input_files <- resolve_manifest_paths(manifest[[file_column]], path)
  source_format <- ifelse(is_vcf_summary_path(input_files), "VCF", "TABLE")
  vcf_sample <- if (is.na(sample_column)) {
    rep(NA_character_, nrow(manifest))
  } else {
    coerce_manifest_text(manifest[[sample_column]])
  }
  if (any(source_format != "VCF" & !is.na(vcf_sample))) {
    stop("Manifest SAMPLE values are only valid for VCF/BCF input files.")
  }
  duplicated_files <- unique(input_files[duplicated(input_files)])
  if (length(duplicated_files) > 0L) {
    for (duplicated_file in duplicated_files) {
      rows <- which(input_files == duplicated_file)
      if (any(source_format[rows] != "VCF") ||
          any(is.na(vcf_sample[rows])) ||
          anyDuplicated(vcf_sample[rows])) {
        stop(
          "A summary-statistic FILE can be repeated only for distinct, explicitly ",
          "named samples in a multi-sample VCF/BCF. Invalid repeated FILE: ",
          duplicated_file
        )
      }
    }
  }

  sample_size <- if (is.na(n_column)) {
    rep(NA_real_, nrow(manifest))
  } else {
    coerce_manifest_number(manifest[[n_column]], "N")
  }
  if (any(!is.na(sample_size) & sample_size <= 0)) {
    stop("Every non-missing LDSC manifest N value must be greater than zero.")
  }
  if (any(source_format == "VCF" & !is.na(sample_size))) {
    stop(
      "Do not provide manifest N for VCF inputs. Binary VCF N is derived as ",
      "NC + NCO; quantitative VCF N is taken directly from NEF."
    )
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
  unmatched_prevalence <- source_format != "VCF" &
    xor(is.na(sample_prev), is.na(population_prev))
  if (any(unmatched_prevalence)) {
    stop(
      "SAMPLE_PREV and POPULATION_PREV must both be NA for quantitative ",
      "traits or both be provided for binary traits. Invalid manifest row(s): ",
      paste(which(unmatched_prevalence) + 1L, collapse = ", "), "."
    )
  }
  vcf_sample_without_population <- source_format == "VCF" &
    !is.na(sample_prev) & is.na(population_prev)
  if (any(vcf_sample_without_population)) {
    stop(
      "A VCF SAMPLE_PREV cannot be supplied without POPULATION_PREV. ",
      "For quantitative VCFs both prevalences must be NA. Invalid manifest ",
      "row(s): ",
      paste(which(vcf_sample_without_population) + 1L, collapse = ", "), "."
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
    source_file = input_files,
    file = input_files,
    source_format = source_format,
    vcf_sample = vcf_sample,
    N = sample_size,
    sample_prev = sample_prev,
    population_prev = population_prev,
    sample_prev_source = ifelse(
      is.na(sample_prev), NA_character_, "manifest"
    ),
    vcf_study_type = NA_character_,
    vcf_has_info = NA
  )
  resolved[["order"]] <- if (is.null(analysis_traits)) {
    seq_len(nrow(resolved))
  } else {
    match(resolved[["trait"]], analysis_traits)
  }
  if (anyNA(resolved[["order"]])) {
    stop("Internal error while reordering the LDSC manifest.")
  }
  data.table::setorder(resolved, order)
  data.table::setcolorder(
    resolved,
    c(
      "order", "trait", "source_file", "file", "source_format",
      "vcf_sample", "N", "sample_prev", "population_prev",
      "sample_prev_source", "vcf_study_type", "vcf_has_info"
    )
  )
  resolved
}

resolve_bcftools_executable <- function(value) {
  requested <- path.expand(trimws(as.character(value)))
  if (length(requested) != 1L || !nzchar(requested)) {
    stop("--bcftools cannot be empty.")
  }
  resolved <- if (grepl("[/\\\\]", requested)) {
    requested
  } else {
    unname(Sys.which(requested))
  }
  if (!nzchar(resolved) || !file.exists(resolved)) {
    stop(
      "VCF input requires bcftools. Install bcftools or provide its executable ",
      "with --bcftools PATH. Requested value: ", requested
    )
  }
  if (file.access(resolved, 1L) != 0L) {
    stop("The bcftools executable is not executable: ", resolved)
  }
  normalizePath(resolved, mustWork = TRUE)
}

bcftools_failure_message <- function(label, status, stderr_file) {
  details <- if (file.exists(stderr_file)) {
    utils::tail(readLines(stderr_file, warn = FALSE), 20L)
  } else {
    character()
  }
  paste0(
    label, " failed with exit status ", status, ".",
    if (length(details) > 0L) {
      paste0("\n  ", paste(details, collapse = "\n  "))
    } else {
      ""
    }
  )
}

run_bcftools_to_file <- function(executable, arguments, output, label) {
  stderr_file <- paste0(output, ".stderr")
  on.exit(unlink(stderr_file, force = TRUE), add = TRUE)
  status <- suppressWarnings(system2(
    executable,
    arguments,
    stdout = output,
    stderr = stderr_file,
    wait = TRUE
  ))
  if (length(status) == 0L || is.null(status)) status <- 0L
  if (!identical(as.integer(status), 0L)) {
    stop(bcftools_failure_message(label, status, stderr_file))
  }
  invisible(output)
}

run_bcftools_append_to_file <- function(executable, arguments, output, label) {
  stderr_file <- paste0(output, ".stderr")
  on.exit(unlink(stderr_file, force = TRUE), add = TRUE)
  command <- paste(
    c(shQuote(executable), arguments, ">>", shQuote(output)),
    collapse = " "
  )
  status <- suppressWarnings(system2(
    "/bin/sh",
    c("-c", shQuote(command)),
    stdout = FALSE,
    stderr = stderr_file,
    wait = TRUE
  ))
  if (length(status) == 0L || is.null(status)) status <- 0L
  if (!identical(as.integer(status), 0L)) {
    stop(bcftools_failure_message(label, status, stderr_file))
  }
  invisible(output)
}

run_bcftools_lines <- function(executable, arguments, label) {
  output <- tempfile("fastasset_bcftools_")
  on.exit(unlink(output, force = TRUE), add = TRUE)
  run_bcftools_to_file(executable, arguments, output, label)
  readLines(output, warn = FALSE)
}

extract_vcf_meta_value <- function(line, key) {
  pattern <- paste0("(^|[,<])", key, "=([^,>]+)")
  matched <- regexec(pattern, line, perl = TRUE)
  values <- regmatches(line, matched)[[1L]]
  if (length(values) < 3L) return(NA_character_)
  sub('^"|"$', "", trimws(values[3L]))
}

vcf_sample_metadata <- function(header, sample) {
  lines <- grep("^##SAMPLE=<", header, value = TRUE)
  if (length(lines) == 0L) return(character())
  ids <- vapply(
    lines, extract_vcf_meta_value, character(1L), key = "ID"
  )
  matching <- lines[!is.na(ids) & ids == sample]
  if (length(matching) == 0L) return(character())
  if (length(matching) > 1L) {
    stop("The VCF has duplicate ##SAMPLE metadata for sample: ", sample)
  }
  matching
}

parse_optional_positive_number <- function(value, label, file) {
  if (length(value) == 0L || is.na(value) || !nzchar(value)) return(NA_real_)
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed) || parsed <= 0) {
    stop("Invalid ", label, " in VCF ##SAMPLE metadata: ", file)
  }
  parsed
}

coerce_vcf_numeric_columns <- function(table, columns, file) {
  for (column in columns) {
    original <- table[[column]]
    parsed <- suppressWarnings(as.numeric(original))
    nonmissing <- !is.na(original) & nzchar(trimws(as.character(original))) &
      trimws(as.character(original)) != "."
    malformed <- nonmissing & !is.finite(parsed)
    if (any(malformed)) {
      stop(
        "VCF FORMAT/", column, " contains malformed numeric value(s) in ",
        file, "."
      )
    }
    table[[column]] <- parsed
  }
  table
}

derive_binary_sample_prevalence <- function(table, supplied_sample_prev,
                                             file) {
  median_cases <- stats::median(table[["NC"]])
  median_controls <- stats::median(table[["NCO"]])
  derived <- median_cases / (median_cases + median_controls)
  source <- paste0(
    "median(FORMAT/NC)/(median(FORMAT/NC)+median(FORMAT/NCO)) ",
    "across valid SNP rows"
  )
  tolerance <- max(1e-6, 0.5 / (median_cases + median_controls))
  if (!is.na(supplied_sample_prev) &&
      abs(supplied_sample_prev - derived) > tolerance) {
    stop(
      "Manifest SAMPLE_PREV does not match the value derived from median ",
      "VCF NC/NCO ",
      "for ", file, ": supplied=", format(supplied_sample_prev, digits = 12),
      ", derived=", format(derived, digits = 12), "."
    )
  }
  list(
    sample_prev = derived,
    median_cases = median_cases,
    median_controls = median_controls,
    source = source
  )
}

standardize_vcf_query_table <- function(table, binary, supplied_sample_prev,
                                        has_info = FALSE,
                                        source_file = "VCF") {
  numeric_columns <- c("POS", "ES", "SE", "LP", "AF")
  if ("NEF" %in% names(table)) numeric_columns <- c(numeric_columns, "NEF")
  if ("SS" %in% names(table)) numeric_columns <- c(numeric_columns, "SS")
  if (binary) numeric_columns <- c(numeric_columns, "NC", "NCO")
  if (has_info) numeric_columns <- c(numeric_columns, "SI")
  table <- coerce_vcf_numeric_columns(
    table, unique(numeric_columns), source_file
  )

  table[["SNP"]] <- trimws(as.character(table[["SNP"]]))
  table[["REF"]] <- toupper(trimws(as.character(table[["REF"]])))
  table[["ALT"]] <- toupper(trimws(as.character(table[["ALT"]])))
  valid <- !is.na(table[["SNP"]]) & nzchar(table[["SNP"]]) &
    table[["SNP"]] != "." &
    grepl("^[ACGT]$", table[["REF"]]) &
    grepl("^[ACGT]$", table[["ALT"]]) &
    table[["REF"]] != table[["ALT"]] &
    is.finite(table[["ES"]]) &
    is.finite(table[["SE"]]) & table[["SE"]] > 0 &
    is.finite(table[["LP"]]) & table[["LP"]] >= 0 &
    is.finite(table[["AF"]]) & table[["AF"]] > 0 & table[["AF"]] < 1
  if (binary) {
    valid <- valid & is.finite(table[["NC"]]) & table[["NC"]] > 0 &
      is.finite(table[["NCO"]]) & table[["NCO"]] > 0
  } else {
    valid <- valid & is.finite(table[["NEF"]]) & table[["NEF"]] > 0
  }
  if (has_info) {
    valid <- valid & is.finite(table[["SI"]]) &
      table[["SI"]] >= 0 & table[["SI"]] <= 1
  }

  input_rows <- nrow(table)
  # Always use explicit row indexing.  `fread()` can be configured globally to
  # return a data.frame; for a data.frame, `table[valid]` means column selection
  # and fails with "undefined columns selected" when `valid` is row-sized.
  table <- table[which(valid), ]
  if (nrow(table) == 0L) {
    stop("No valid biallelic SNP rows remained after VCF conversion: ", source_file)
  }
  if (anyDuplicated(table[["SNP"]])) {
    duplicates <- unique(table[["SNP"]][duplicated(table[["SNP"]])])
    stop(
      "VCF contains duplicate SNP identifiers after filtering: ",
      paste(utils::head(duplicates, 10L), collapse = ", "),
      if (length(duplicates) > 10L) " ..." else "",
      ". Normalize/de-duplicate the VCF before analysis: ", source_file
    )
  }

  minimum_log10 <- -log10(.Machine$double.xmin)
  p_was_capped <- table[["LP"]] > minimum_log10
  raw_p <- 10^(-pmin(table[["LP"]], minimum_log10))

  prevalence <- if (binary) {
    derive_binary_sample_prevalence(
      table, supplied_sample_prev, source_file
    )
  } else {
    list(
      sample_prev = NA_real_, median_cases = NA_real_,
      median_controls = NA_real_,
      source = "not applicable (quantitative)"
    )
  }
  per_snp_n <- if (binary) {
    table[["NC"]] + table[["NCO"]]
  } else {
    table[["NEF"]]
  }
  if (binary && "SS" %in% names(table)) {
    comparable <- is.finite(table[["SS"]])
    inconsistent <- comparable & abs(table[["SS"]] - per_snp_n) > 0.5
    if (any(inconsistent)) {
      stop(
        "VCF FORMAT/SS is inconsistent with FORMAT/NC + FORMAT/NCO for ",
        sum(inconsistent), " converted row(s): ", source_file
      )
    }
  }

  output <- data.table::data.table(
    SNP = table[["SNP"]],
    A1 = table[["ALT"]],
    A2 = table[["REF"]],
    BETA = table[["ES"]],
    SE = table[["SE"]],
    P = raw_p,
    N = per_snp_n,
    MAF = pmin(table[["AF"]], 1 - table[["AF"]])
  )
  if (has_info) output[["INFO"]] <- table[["SI"]]
  output[["CHR"]] <- as.character(table[["CHR"]])
  output[["POS"]] <- table[["POS"]]
  output[["LP"]] <- table[["LP"]]
  output[["AF"]] <- table[["AF"]]
  if ("SS" %in% names(table)) output[["SS"]] <- table[["SS"]]
  if (binary) {
    output[["NC"]] <- table[["NC"]]
    output[["NCO"]] <- table[["NCO"]]
  }
  if ("NEF" %in% names(table)) output[["NEF"]] <- table[["NEF"]]

  list(
    table = output,
    sample_prev = prevalence$sample_prev,
    sample_prev_source = prevalence$source,
    median_cases = prevalence$median_cases,
    median_controls = prevalence$median_controls,
    input_rows = input_rows,
    output_rows = nrow(output),
    dropped_rows = input_rows - nrow(output),
    p_values_capped = sum(p_was_capped)
  )
}

write_generated_table <- function(table, path) {
  temporary <- paste0(path, ".building-", Sys.getpid(), ".gz")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  data.table::fwrite(
    table, temporary, sep = "\t", quote = FALSE, na = "NA",
    compress = "gzip"
  )
  if (!file.rename(temporary, path)) {
    stop("Could not move generated VCF-conversion table into place: ", path)
  }
  invisible(path)
}

write_conversion_metadata <- function(values, path) {
  metadata <- data.table::data.table(
    parameter = names(values),
    value = vapply(values, as.character, character(1L))
  )
  temporary <- paste0(path, ".building-", Sys.getpid())
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  data.table::fwrite(metadata, temporary, sep = "\t", quote = FALSE)
  if (!file.rename(temporary, path)) {
    stop("Could not move VCF-conversion metadata into place: ", path)
  }
  invisible(path)
}

read_conversion_metadata <- function(path) {
  metadata <- data.table::fread(path, colClasses = "character")
  if (!all(c("parameter", "value") %in% names(metadata)) ||
      anyDuplicated(metadata[["parameter"]])) {
    stop("Malformed VCF-conversion metadata: ", path)
  }
  stats::setNames(as.list(metadata[["value"]]), metadata[["parameter"]])
}

read_bcftools_query_output <- function(path, expected_columns, source_file) {
  query <- data.table::fread(
    path,
    header = TRUE,
    na.strings = c(".", "NA", ""),
    check.names = FALSE,
    data.table = TRUE
  )
  if (ncol(query) != length(expected_columns)) {
    stop(
      "bcftools query produced ", ncol(query), " column(s), but ",
      length(expected_columns), " were expected for ", source_file,
      ". Expected columns: ", paste(expected_columns, collapse = ", "), "."
    )
  }
  if (!identical(names(query), expected_columns)) {
    stop(
      "Internal bcftools-query header mismatch for ", source_file,
      ". Observed: ", paste(names(query), collapse = ", "),
      ". Expected: ", paste(expected_columns, collapse = ", "), "."
    )
  }
  query
}

prepare_one_vcf_summary <- function(source_file, requested_sample,
                                    population_prev, supplied_sample_prev,
                                    output_file, metadata_file, bcftools,
                                    bcftools_version) {
  required_metadata <- c(
    "sample", "study_type", "sample_prev", "sample_prev_source",
    "prevalence_method", "has_info"
  )
  if (file.exists(output_file) && file.exists(metadata_file)) {
    metadata <- read_conversion_metadata(metadata_file)
    if (all(required_metadata %in% names(metadata))) {
      return(list(
        file = normalizePath(output_file, mustWork = TRUE),
        sample = metadata[["sample"]],
        study_type = metadata[["study_type"]],
        sample_prev = suppressWarnings(as.numeric(metadata[["sample_prev"]])),
        sample_prev_source = metadata[["sample_prev_source"]],
        has_info = identical(metadata[["has_info"]], "TRUE")
      ))
    }
  }

  header <- run_bcftools_lines(
    bcftools,
    c("view", "--header-only", shQuote(source_file)),
    paste0("bcftools header inspection for ", source_file)
  )
  samples <- run_bcftools_lines(
    bcftools,
    c("query", "--list-samples", shQuote(source_file)),
    paste0("bcftools sample inspection for ", source_file)
  )
  samples <- samples[nzchar(samples)]
  if (length(samples) == 0L) {
    stop("VCF has no sample column: ", source_file)
  }
  sample <- requested_sample
  if (is.na(sample)) {
    if (length(samples) != 1L) {
      stop(
        "VCF contains multiple samples; provide the manifest SAMPLE value for ",
        source_file, ". Available samples: ", paste(samples, collapse = ", ")
      )
    }
    sample <- samples[1L]
  } else if (!sample %in% samples) {
    stop(
      "Manifest SAMPLE is not present in VCF ", source_file, ": ", sample,
      ". Available samples: ", paste(samples, collapse = ", ")
    )
  }

  sample_metadata <- vcf_sample_metadata(header, sample)
  study_type <- if (length(sample_metadata) == 1L) {
    extract_vcf_meta_value(sample_metadata, "StudyType")
  } else {
    NA_character_
  }
  normalized_study_type <- toupper(gsub("[^A-Za-z]", "", study_type))
  header_binary <- normalized_study_type %in% c("CASECONTROL", "BINARY")
  header_quantitative <- normalized_study_type %in%
    c("CONTINUOUS", "QUANTITATIVE")
  binary <- !is.na(population_prev)
  if (header_binary && !binary) {
    stop(
      "VCF ##SAMPLE StudyType is CaseControl, but POPULATION_PREV is missing ",
      "from the manifest: ", source_file
    )
  }
  if (header_quantitative && binary) {
    stop(
      "VCF ##SAMPLE StudyType is Continuous, but POPULATION_PREV was supplied: ",
      source_file
    )
  }
  if (is.na(study_type)) {
    study_type <- if (binary) "CaseControl (inferred)" else "Continuous (inferred)"
  }

  format_lines <- grep("^##FORMAT=<ID=", header, value = TRUE)
  format_ids <- sub("^##FORMAT=<ID=([^,>]+).*$", "\\1", format_lines)
  required_fields <- c("ES", "SE", "LP", "AF")
  if (binary) {
    required_fields <- c(required_fields, "NC", "NCO")
  } else {
    required_fields <- c(required_fields, "NEF")
  }
  missing_fields <- setdiff(required_fields, format_ids)
  if (length(missing_fields) > 0L) {
    stop(
      "VCF is missing required FORMAT field(s): ",
      paste(missing_fields, collapse = ", "), ". File: ", source_file
    )
  }
  has_info <- "SI" %in% format_ids
  optional_fields <- intersect(c("SS", "SI", "NEF"), format_ids)
  queried_fields <- c(required_fields, optional_fields)
  queried_fields <- queried_fields[!duplicated(queried_fields)]
  query_format <- paste0(
    paste(c("%CHROM", "%POS", "%ID", "%REF", "%ALT"), collapse = "\\t"),
    "[\\t", paste0("%", queried_fields, collapse = "\\t"), "]\\n"
  )
  raw_output <- paste0(output_file, ".bcftools-query-", Sys.getpid(), ".tsv")
  on.exit(unlink(raw_output, force = TRUE), add = TRUE)
  expected_columns <- c("CHR", "POS", "SNP", "REF", "ALT", queried_fields)
  writeLines(
    paste(expected_columns, collapse = "\t"),
    raw_output,
    useBytes = TRUE
  )
  run_bcftools_append_to_file(
    bcftools,
    c(
      "query", "--samples", shQuote(sample), "--format",
      shQuote(query_format), shQuote(source_file)
    ),
    raw_output,
    paste0("bcftools query for ", source_file)
  )
  query <- read_bcftools_query_output(
    raw_output, expected_columns, source_file
  )

  header_cases <- if (length(sample_metadata) == 1L) {
    parse_optional_positive_number(
      extract_vcf_meta_value(sample_metadata, "TotalCases"),
      "TotalCases", source_file
    )
  } else {
    NA_real_
  }
  header_controls <- if (length(sample_metadata) == 1L) {
    parse_optional_positive_number(
      extract_vcf_meta_value(sample_metadata, "TotalControls"),
      "TotalControls", source_file
    )
  } else {
    NA_real_
  }
  if (xor(is.finite(header_cases), is.finite(header_controls))) {
    stop(
      "VCF ##SAMPLE metadata must provide both TotalCases and TotalControls ",
      "or neither: ", source_file
    )
  }

  standardized <- standardize_vcf_query_table(
    query,
    binary = binary,
    supplied_sample_prev = supplied_sample_prev,
    has_info = has_info,
    source_file = source_file
  )
  write_generated_table(standardized$table, output_file)
  write_conversion_metadata(
    list(
      source_file = source_file,
      converted_file = output_file,
      sample = sample,
      study_type = study_type,
      trait_type = if (binary) "binary" else "quantitative",
      population_prev = if (binary) population_prev else NA_real_,
      sample_prev = standardized$sample_prev,
      sample_prev_source = standardized$sample_prev_source,
      prevalence_method = if (binary) {
        "ratio of separate per-SNP NC and NCO medians"
      } else {
        "not applicable"
      },
      median_cases = standardized$median_cases,
      median_controls = standardized$median_controls,
      header_total_cases = header_cases,
      header_total_controls = header_controls,
      has_info = has_info,
      input_rows = standardized$input_rows,
      output_rows = standardized$output_rows,
      dropped_rows = standardized$dropped_rows,
      p_values_capped = standardized$p_values_capped,
      p_conversion = "P=10^(-LP); values below .Machine$double.xmin capped",
      effect_allele = "A1=ALT; BETA=FORMAT/ES",
      other_allele = "A2=REF",
      ldsc_n = if (binary) "FORMAT/NC + FORMAT/NCO" else "FORMAT/NEF",
      bcftools = bcftools,
      bcftools_version = bcftools_version
    ),
    metadata_file
  )
  list(
    file = normalizePath(output_file, mustWork = TRUE),
    sample = sample,
    study_type = study_type,
    sample_prev = standardized$sample_prev,
    sample_prev_source = standardized$sample_prev_source,
    has_info = has_info
  )
}

prepare_manifest_vcf_files <- function(params, manifest, generation_dir) {
  rows <- which(manifest[["source_format"]] == "VCF")
  if (length(rows) == 0L) {
    return(list(
      manifest = manifest,
      bcftools = NA_character_,
      bcftools_version = NA_character_
    ))
  }
  bcftools <- resolve_bcftools_executable(params$bcftools)
  version_lines <- run_bcftools_lines(
    bcftools, "--version", "bcftools version check"
  )
  bcftools_version <- if (length(version_lines) > 0L) {
    version_lines[1L]
  } else {
    "unknown"
  }
  conversion_dir <- file.path(generation_dir, "vcf_converted")
  dir.create(conversion_dir, recursive = TRUE, showWarnings = FALSE)
  output_files <- file.path(
    conversion_dir, sprintf("trait_%04d.tsv.gz", manifest[["order"]][rows])
  )
  metadata_files <- file.path(
    conversion_dir, sprintf("trait_%04d.metadata.tsv", manifest[["order"]][rows])
  )
  worker <- function(task) {
    row <- rows[task]
    tryCatch(
      list(
        ok = TRUE,
        result = prepare_one_vcf_summary(
          source_file = manifest[["source_file"]][row],
          requested_sample = manifest[["vcf_sample"]][row],
          population_prev = manifest[["population_prev"]][row],
          supplied_sample_prev = manifest[["sample_prev"]][row],
          output_file = output_files[task],
          metadata_file = metadata_files[task],
          bcftools = bcftools,
          bcftools_version = bcftools_version
        )
      ),
      error = function(error) list(ok = FALSE, error = conditionMessage(error))
    )
  }
  workers <- min(params$vcf_cores, length(rows))
  message(
    "Converting ", length(rows), " VCF summary-statistic file(s) with ",
    workers, " bcftools worker(s)."
  )
  results <- if (.Platform$OS.type == "windows" || workers == 1L) {
    lapply(seq_along(rows), worker)
  } else {
    parallel::mclapply(
      seq_along(rows), worker, mc.cores = workers, mc.preschedule = FALSE
    )
  }
  failed <- which(!vapply(
    results,
    function(value) is.list(value) && isTRUE(value$ok),
    logical(1L)
  ))
  if (length(failed) > 0L) {
    messages <- vapply(
      failed,
      function(index) {
        detail <- if (is.list(results[[index]]) &&
                      !is.null(results[[index]]$error)) {
          results[[index]]$error
        } else {
          "worker terminated without returning an error message"
        }
        paste0(manifest[["trait"]][rows[index]], ": ", detail)
      },
      character(1L)
    )
    stop("VCF conversion failed:\n  ", paste(messages, collapse = "\n  "))
  }
  for (task in seq_along(rows)) {
    row <- rows[task]
    result <- results[[task]]$result
    data.table::set(manifest, i = row, j = "file", value = result$file)
    data.table::set(manifest, i = row, j = "vcf_sample", value = result$sample)
    data.table::set(
      manifest, i = row, j = "vcf_study_type", value = result$study_type
    )
    data.table::set(
      manifest, i = row, j = "sample_prev", value = result$sample_prev
    )
    data.table::set(
      manifest, i = row, j = "sample_prev_source",
      value = result$sample_prev_source
    )
    data.table::set(
      manifest, i = row, j = "vcf_has_info", value = result$has_info
    )
  }
  list(
    manifest = manifest,
    bcftools = bcftools,
    bcftools_version = bcftools_version
  )
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
    "ldsc_generation_version=2026-08-15-v3-median-prevalence",
    paste0("traits=", paste(manifest$trait, collapse = ";")),
    paste0("source_format=", paste(manifest$source_format, collapse = ";")),
    paste0("vcf_sample=", paste(manifest$vcf_sample, collapse = ";")),
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
    paste0("bcftools=", params$bcftools),
    paste0("vcf_cores=", params$vcf_cores),
    file_metadata_lines("manifest", params$sumstats_manifest),
    file_metadata_lines("hm3", params$hm3),
    file_metadata_lines("summary", manifest$source_file),
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

generate_ldsc_object <- function(params, analysis_traits, prepared_input) {
  if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
    stop(
      "Automatic-LDSC mode requires the GenomicSEM package. Re-run install.sh ",
      "without FASTASSET_SKIP_DEPENDENCIES=1 or install GenomicSEM first."
    )
  }

  manifest <- data.table::copy(prepared_input$manifest)
  if (!identical(as.character(manifest[["trait"]]), as.character(analysis_traits))) {
    stop(
      "Internal trait-order mismatch between the prepared manifest input and ",
      "automatic LDSC generation."
    )
  }
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

  needs_munge <- !file.exists(manifest$munged_file)
  run_munge_group <- function(rows, column_map, label) {
    rows <- rows[needs_munge[rows]]
    if (length(rows) == 0L) return(invisible(NULL))
    workers <- min(params$munge_cores, length(rows))
    message(
      "Munging ", length(rows), " ", label,
      " summary-statistic file(s) with ", workers, " worker(s)."
    )
    GenomicSEM::munge(
      files = manifest$file[rows],
      hm3 = params$hm3,
      trait.names = manifest$munged_prefix[rows],
      N = manifest$N[rows],
      info.filter = params$munge_info_filter,
      maf.filter = params$munge_maf_filter,
      log.name = file.path(
        generation_dir, paste0(params$run_name, "_munge_", label)
      ),
      column.names = column_map,
      parallel = length(rows) > 1L && params$munge_cores > 1L,
      cores = workers,
      overwrite = FALSE
    )
    invisible(NULL)
  }
  if (any(needs_munge)) {
    table_rows <- which(manifest$source_format == "TABLE")
    vcf_info_rows <- which(
      manifest$source_format == "VCF" & manifest$vcf_has_info
    )
    vcf_no_info_rows <- which(
      manifest$source_format == "VCF" & !manifest$vcf_has_info
    )
    vcf_map <- list(
      SNP = "SNP", A1 = "A1", A2 = "A2", effect = "BETA",
      P = "P", N = "N", MAF = "MAF"
    )
    genomicsem_map <- params$munge_column_map[
      names(params$munge_column_map) %in%
        c("SNP", "A1", "A2", "effect", "INFO", "P", "N", "MAF", "Z")
    ]
    run_munge_group(table_rows, genomicsem_map, "table")
    run_munge_group(
      vcf_info_rows, c(vcf_map, list(INFO = "INFO")), "vcf_info"
    )
    run_munge_group(vcf_no_info_rows, vcf_map, "vcf_no_info")
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
      "bcftools", "bcftools_version", "vcf_cores", "number_vcf_inputs",
      "ldsc_chisq_max", "ldsc_object_name", "generated_rdata"
    ),
    value = as.character(c(
      signature, number_traits, params$ldsc_blocks, effective_blocks,
      params$ldsc_chr, params$ld_ref, params$wld_ref, params$hm3,
      params$sumstats_manifest, params$munge_info_filter,
      params$munge_maf_filter, params$munge_cores,
      prepared_input$bcftools, prepared_input$bcftools_version,
      params$vcf_cores, sum(manifest$source_format == "VCF"),
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

resolve_ldsc_source <- function(params, analysis_traits, prepared_input) {
  if (params$ldsc_mode == "existing") {
    return(list(
      path = params$ldsc_rdata,
      generation_dir = NA_character_,
      signature = NA_character_,
      reused = NA
    ))
  }
  generate_ldsc_object(params, analysis_traits, prepared_input)
}
