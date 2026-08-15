summary_column_from_map <- function(column_names, mapping, mapping_key,
                                    aliases, label, required = TRUE) {
  mapped <- mapping[[mapping_key]]
  if (!is.null(mapped)) {
    if (!mapped %in% column_names) {
      stop(
        "--munge-column-map assigns ", label, " to '", mapped,
        "', but that column is absent. Available columns: ",
        paste(column_names, collapse = ", "), "."
      )
    }
    return(mapped)
  }

  normalized <- normalize_manifest_header(column_names)
  if (anyDuplicated(normalized)) {
    stop(
      "Summary-statistic table has column names that become duplicates after ",
      "normalization: ", paste(column_names[duplicated(normalized)], collapse = ", "),
      "."
    )
  }
  for (alias in aliases) {
    hit <- which(normalized == alias)
    if (length(hit) == 1L) return(column_names[hit])
  }
  if (!required) return(NA_character_)
  stop(
    "Could not identify the ", label, " column. Accepted automatic names: ",
    paste(aliases, collapse = ", "),
    ". Supply an explicit --munge-column-map entry."
  )
}

coerce_summary_numeric <- function(values, label, file) {
  text <- trimws(as.character(values))
  missing <- is.na(values) | !nzchar(text) |
    toupper(text) %in% c("NA", "N/A", "NULL", ".")
  parsed <- suppressWarnings(as.numeric(text))
  malformed <- !missing & !is.finite(parsed)
  if (any(malformed)) {
    stop(
      label, " contains malformed or non-finite value(s) in ", file,
      ". First affected input row(s): ",
      paste(utils::head(which(malformed) + 1L, 10L), collapse = ", "), "."
    )
  }
  parsed[missing] <- NA_real_
  parsed
}

read_manifest_trait_for_fastasset <- function(manifest_row, params) {
  file <- as.character(manifest_row[["file"]][1L])
  trait <- as.character(manifest_row[["trait"]][1L])
  if (length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("Internal error: the resolved manifest row has no prepared FILE value.")
  }
  binary <- !is.na(manifest_row[["population_prev"]][1L])
  source_format <- manifest_row[["source_format"]][1L]
  header <- names(data.table::fread(
    file, nrows = 0L, check.names = FALSE, showProgress = FALSE
  ))
  if (length(header) == 0L) stop("Summary-statistic file has no header: ", file)

  mapping <- if (identical(source_format, "VCF")) list() else {
    params$munge_column_map
  }
  snp_column <- summary_column_from_map(
    header, mapping, "SNP",
    c("SNP", "ID", "RSID", "MARKERNAME", "MARKER_NAME"), "SNP/ID"
  )
  a1_column <- summary_column_from_map(
    header, mapping, "A1",
    c("A1", "EA", "EFFECT_ALLELE", "ALT"), "effect allele (A1)"
  )
  a2_column <- summary_column_from_map(
    header, mapping, "A2",
    c("A2", "NEA", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "REF"),
    "other allele (A2)"
  )
  beta_column <- summary_column_from_map(
    header, mapping, "effect",
    c("BETA", "ES", "EFFECT", "B", "LOG_OR", "LOGOR"), "effect/BETA"
  )
  se_column <- summary_column_from_map(
    header, mapping, "SE",
    c("SE", "STANDARD_ERROR", "STDERR", "SEBETA", "SE_BETA"), "SE"
  )

  ncase_column <- if (binary) {
    summary_column_from_map(
      header, mapping, "NCASE", c("NC", "N_CASE", "NCASE", "CASES"),
      "case count", required = FALSE
    )
  } else {
    NA_character_
  }
  ncontrol_column <- if (binary) {
    summary_column_from_map(
      header, mapping, "NCONTROL",
      c("NCO", "N_CONTROL", "NCONTROL", "CONTROLS"),
      "control count", required = FALSE
    )
  } else {
    NA_character_
  }
  if (binary && xor(is.na(ncase_column), is.na(ncontrol_column))) {
    stop(
      "Binary sample-size construction requires both case and control columns, ",
      "not only one. Trait: ", trait, "; file: ", file, "."
    )
  }

  n_column <- summary_column_from_map(
    header, mapping, "N",
    if (binary) {
      c("N", "SS", "TOTAL_N", "SAMPLE_SIZE")
    } else {
      c("NEF", "N", "N_EFF", "NEFF", "SS", "TOTAL_N", "SAMPLE_SIZE")
    },
    if (binary) "total sample size" else "quantitative NEF/total N",
    required = FALSE
  )
  manifest_n <- manifest_row[["N"]][1L]
  if (!binary && is.na(n_column) && is.na(manifest_n)) {
    stop(
      "Quantitative trait '", trait, "' has no per-SNP NEF/N column and no ",
      "constant manifest N: ", file, "."
    )
  }
  if (binary && is.na(ncase_column) && is.na(n_column) && is.na(manifest_n)) {
    stop(
      "Binary trait '", trait, "' requires NC and NCO columns, or a total-N ",
      "column, or a constant manifest N: ", file, "."
    )
  }

  selected <- c(
    snp_column, a1_column, a2_column, beta_column, se_column,
    ncase_column, ncontrol_column, n_column
  )
  selected <- unique(selected[!is.na(selected)])
  table <- data.table::fread(
    file,
    select = selected,
    na.strings = c("", "NA", "N/A", "NULL", "."),
    check.names = FALSE,
    showProgress = FALSE
  )
  id <- trimws(as.character(table[[snp_column]]))
  a1 <- toupper(trimws(as.character(table[[a1_column]])))
  a2 <- toupper(trimws(as.character(table[[a2_column]])))
  invalid_identifier <- is.na(id) | !nzchar(id) | id == "."
  invalid_allele <- is.na(a1) | is.na(a2) |
    !grepl("^[ACGT]$", a1) | !grepl("^[ACGT]$", a2) | a1 == a2
  if (any(invalid_identifier) || any(invalid_allele)) {
    stop(
      "FastASSET input construction requires non-missing SNP IDs and distinct ",
      "single-base A/C/G/T alleles. Trait '", trait, "' has ",
      sum(invalid_identifier), " invalid ID row(s) and ",
      sum(invalid_allele), " invalid allele row(s): ", file, "."
    )
  }
  if (anyDuplicated(id)) {
    duplicates <- unique(id[duplicated(id)])
    stop(
      "Trait '", trait, "' contains duplicate SNP IDs: ",
      paste(utils::head(duplicates, 10L), collapse = ", "),
      if (length(duplicates) > 10L) " ..." else "", ". File: ", file
    )
  }

  beta <- coerce_summary_numeric(table[[beta_column]], "BETA", file)
  se <- coerce_summary_numeric(table[[se_column]], "SE", file)
  if (!is.na(ncase_column)) {
    ncase <- coerce_summary_numeric(table[[ncase_column]], "NCASE", file)
    ncontrol <- coerce_summary_numeric(
      table[[ncontrol_column]], "NCONTROL", file
    )
    neff <- ncase * ncontrol / (ncase + ncontrol)
    sample_size_source <- paste0(
      ncase_column, "*", ncontrol_column, "/(", ncase_column, "+",
      ncontrol_column, ")"
    )
  } else {
    total_n <- if (!is.na(n_column)) {
      coerce_summary_numeric(table[[n_column]], "N", file)
    } else {
      rep(manifest_n, nrow(table))
    }
    if (binary) {
      sample_prev <- manifest_row[["sample_prev"]][1L]
      if (!is.finite(sample_prev)) {
        stop(
          "Binary trait '", trait,
          "' needs SAMPLE_PREV when NC/NCO are unavailable: ", file, "."
        )
      }
      neff <- total_n * sample_prev * (1 - sample_prev)
      sample_size_source <- paste0(
        if (!is.na(n_column)) n_column else "manifest N",
        "*SAMPLE_PREV*(1-SAMPLE_PREV)"
      )
    } else {
      neff <- total_n
      sample_size_source <- if (!is.na(n_column)) n_column else "manifest N"
    }
  }

  list(
    table = data.table::data.table(
      ID = id, A1 = a1, A2 = a2,
      Beta = beta, SE = se, NEF = neff
    ),
    audit = data.table::data.table(
      order = manifest_row[["order"]][1L],
      trait = trait,
      source_file = manifest_row[["source_file"]][1L],
      prepared_file = file,
      source_format = source_format,
      trait_type = if (binary) "binary" else "quantitative",
      input_rows = nrow(table),
      snp_column = snp_column,
      a1_column = a1_column,
      a2_column = a2_column,
      beta_column = beta_column,
      se_column = se_column,
      fastasset_neff_source = sample_size_source
    )
  )
}

append_fastasset_trait <- function(master, incoming, trait) {
  exact_matches <- 0L
  swapped_matches <- 0L
  incompatible_matches <- 0L
  incompatible_positions <- integer()
  failed_alignments <- data.table::data.table(
    trait = character(),
    ID = character(),
    input_row = integer(),
    reference_A1 = character(),
    reference_A2 = character(),
    incoming_A1 = character(),
    incoming_A2 = character(),
    reason = character(),
    action = character()
  )
  if (is.null(master)) {
    master <- data.table::data.table(
      ID = incoming[["ID"]],
      reference_A1 = incoming[["A1"]],
      reference_A2 = incoming[["A2"]]
    )
    index <- seq_len(nrow(incoming))
    new_snps <- nrow(incoming)
    flipped <- rep(FALSE, nrow(incoming))
  } else {
    index <- match(incoming[["ID"]], master[["ID"]])
    existing <- !is.na(index)
    flipped <- rep(FALSE, nrow(incoming))
    if (any(existing)) {
      existing_positions <- which(existing)
      reference_a1 <- master[["reference_A1"]][index[existing]]
      reference_a2 <- master[["reference_A2"]][index[existing]]
      exact <- incoming[["A1"]][existing] == reference_a1 &
        incoming[["A2"]][existing] == reference_a2
      swapped <- incoming[["A1"]][existing] == reference_a2 &
        incoming[["A2"]][existing] == reference_a1
      incompatible <- !(exact | swapped)
      exact_matches <- sum(exact)
      swapped_matches <- sum(swapped)
      incompatible_matches <- sum(incompatible)
      if (any(incompatible)) {
        incompatible_positions <- existing_positions[incompatible]
        failed_alignments <- data.table::data.table(
          trait = rep(as.character(trait), incompatible_matches),
          ID = incoming[["ID"]][incompatible_positions],
          input_row = as.integer(incompatible_positions),
          reference_A1 = reference_a1[incompatible],
          reference_A2 = reference_a2[incompatible],
          incoming_A1 = incoming[["A1"]][incompatible_positions],
          incoming_A2 = incoming[["A2"]][incompatible_positions],
          reason = rep("INCOMPATIBLE_ALLELES", incompatible_matches),
          action = rep(
            "BETA_SE_NEF_SET_TO_NA_TRAIT_EXCLUDED_FOR_SNP",
            incompatible_matches
          )
        )
      }
      flipped[existing_positions[swapped]] <- TRUE
    }

    novel <- which(is.na(index))
    new_snps <- length(novel)
    if (new_snps > 0L) {
      additions <- data.table::data.table(
        ID = incoming[["ID"]][novel],
        reference_A1 = incoming[["A1"]][novel],
        reference_A2 = incoming[["A2"]][novel]
      )
      master <- data.table::rbindlist(
        list(master, additions), use.names = TRUE, fill = TRUE
      )
      index[novel] <- nrow(master) - new_snps + seq_len(new_snps)
    }
  }

  beta <- incoming[["Beta"]]
  standard_error <- incoming[["SE"]]
  neff <- incoming[["NEF"]]
  beta[flipped] <- -beta[flipped]
  if (length(incompatible_positions) > 0L) {
    beta[incompatible_positions] <- NA_real_
    standard_error[incompatible_positions] <- NA_real_
    neff[incompatible_positions] <- NA_real_
  }
  data.table::set(master, i = index, j = paste0(trait, ".Beta"), value = beta)
  data.table::set(
    master, i = index, j = paste0(trait, ".SE"), value = standard_error
  )
  data.table::set(master, i = index, j = paste0(trait, ".NEF"), value = neff)
  list(
    master = master,
    allele_flips = sum(flipped),
    new_snps = new_snps,
    exact_matches = exact_matches,
    swapped_matches = swapped_matches,
    incompatible_matches = incompatible_matches,
    failed_alignments = failed_alignments
  )
}

build_manifest_fastasset_input <- function(manifest, params, output_file,
                                           audit_file, allele_file,
                                           failed_file) {
  master <- NULL
  audit_rows <- vector("list", nrow(manifest))
  failed_rows <- vector("list", nrow(manifest))
  for (row in seq_len(nrow(manifest))) {
    manifest_row <- lapply(manifest, function(column) column[row])
    trait_input <- read_manifest_trait_for_fastasset(manifest_row, params)
    appended <- append_fastasset_trait(
      master, trait_input$table, manifest[["trait"]][row]
    )
    master <- appended$master
    trait_input$audit[["allele_flips_to_reference"]] <- appended$allele_flips
    trait_input$audit[["new_union_snps"]] <- appended$new_snps
    trait_input$audit[["reference_established_snps"]] <- appended$new_snps
    trait_input$audit[["exact_allele_matches_to_reference"]] <-
      appended$exact_matches
    trait_input$audit[["swapped_allele_matches_to_reference"]] <-
      appended$swapped_matches
    trait_input$audit[["incompatible_allele_matches"]] <-
      appended$incompatible_matches
    audit_rows[[row]] <- trait_input$audit
    failed_rows[[row]] <- data.table::copy(appended$failed_alignments)
    failed_rows[[row]][["source_file"]] <- rep(
      as.character(manifest[["source_file"]][row]),
      nrow(failed_rows[[row]])
    )
    data.table::setcolorder(
      failed_rows[[row]],
      c(
        "trait", "ID", "source_file", "input_row",
        "reference_A1", "reference_A2", "incoming_A1", "incoming_A2",
        "reason", "action"
      )
    )
    if (row == 1L || row == nrow(manifest) || row %% 10L == 0L) {
      message(
        "Prepared FastASSET columns for trait ", row, "/", nrow(manifest),
        ": ", manifest[["trait"]][row]
      )
    }
  }
  if (is.null(master) || nrow(master) == 0L) {
    stop("No SNP rows were available for FastASSET input construction.")
  }

  final_union_snps <- nrow(master)
  for (row in seq_along(audit_rows)) {
    audit_rows[[row]][["final_union_snps"]] <- final_union_snps
    audit_rows[[row]][["snps_absent_from_trait"]] <-
      final_union_snps - audit_rows[[row]][["input_rows"]]
  }

  data.table::setorder(master, ID)
  allele_reference <- master[, c("ID", "reference_A1", "reference_A2"), with = FALSE]
  analysis_columns <- c(
    "ID",
    unlist(lapply(manifest[["trait"]], function(trait) {
      paste0(trait, c(".Beta", ".SE", ".NEF"))
    }), use.names = FALSE)
  )
  analysis_input <- master[, analysis_columns, with = FALSE]
  write_generated_table(analysis_input, output_file)
  write_generated_table(allele_reference, allele_file)
  data.table::fwrite(
    data.table::rbindlist(audit_rows, use.names = TRUE, fill = TRUE),
    audit_file, sep = "\t", quote = FALSE, na = "NA"
  )
  failed_table <- data.table::rbindlist(
    failed_rows, use.names = TRUE, fill = TRUE
  )
  data.table::fwrite(
    failed_table, failed_file, sep = "\t", quote = FALSE, na = "NA"
  )
  if (nrow(failed_table) > 0L) {
    message(
      "Excluded ", nrow(failed_table),
      " incompatible SNP-trait observation(s); details: ", failed_file
    )
  }
  list(
    path = normalizePath(output_file, mustWork = TRUE),
    n_snps = nrow(master),
    n_failed_alignments = nrow(failed_table),
    failed_alignment_file = normalizePath(failed_file, mustWork = TRUE)
  )
}

manifest_input_signature <- function(params, manifest) {
  mapping <- if (length(params$munge_column_map) == 0L) {
    "AUTO"
  } else {
    paste(
      names(params$munge_column_map),
      unlist(params$munge_column_map, use.names = FALSE),
      sep = "=", collapse = ","
    )
  }
  make_text_signature(c(
    "manifest_input_version=2026-08-15-v4-recoverable-allele-failures",
    paste0("traits=", paste(manifest[["trait"]], collapse = ";")),
    paste0("source_format=", paste(manifest[["source_format"]], collapse = ";")),
    paste0("vcf_sample=", paste(manifest[["vcf_sample"]], collapse = ";")),
    paste0("N=", paste(manifest[["N"]], collapse = ";")),
    paste0("sample_prev=", paste(manifest[["sample_prev"]], collapse = ";")),
    paste0(
      "population_prev=", paste(manifest[["population_prev"]], collapse = ";")
    ),
    paste0("column_map=", mapping),
    paste0("bcftools=", params$bcftools),
    file_metadata_lines("manifest", params$sumstats_manifest),
    file_metadata_lines("summary", manifest[["source_file"]])
  ))
}

read_preparation_provenance <- function(path) {
  table <- data.table::fread(path, colClasses = "character")
  if (!all(c("parameter", "value") %in% names(table))) {
    stop("Malformed input-preparation provenance: ", path)
  }
  stats::setNames(as.list(table[["value"]]), table[["parameter"]])
}

prepare_manifest_analysis_input <- function(params) {
  manifest <- read_ldsc_manifest(params$sumstats_manifest)
  signature <- manifest_input_signature(params, manifest)
  preparation_dir <- file.path(
    params$output_dir, paste0(params$run_name, "_prepared_input_", signature)
  )
  dir.create(preparation_dir, recursive = TRUE, showWarnings = FALSE)
  analysis_file <- file.path(
    preparation_dir, paste0(params$run_name, "_fastasset_wide.tsv.gz")
  )
  allele_file <- file.path(preparation_dir, "fastasset_allele_reference.tsv.gz")
  audit_file <- file.path(preparation_dir, "fastasset_input_build_audit.tsv")
  failed_file <- file.path(
    preparation_dir, "fastasset_failed_allele_alignments.tsv"
  )
  resolved_manifest_file <- file.path(preparation_dir, "manifest_resolved.tsv")
  provenance_file <- file.path(preparation_dir, "input_preparation_provenance.tsv")
  completion_file <- file.path(preparation_dir, "input_preparation.complete")

  complete <- all(file.exists(c(
    analysis_file, allele_file, audit_file, failed_file, resolved_manifest_file,
    provenance_file, completion_file
  ))) && identical(
    trimws(readLines(completion_file, warn = FALSE)[1L]), signature
  )
  if (complete) {
    provenance <- read_preparation_provenance(provenance_file)
    message("Validated manifest-derived FastASSET input exists; reusing: ", analysis_file)
    return(list(
      path = normalizePath(analysis_file, mustWork = TRUE),
      manifest = data.table::fread(
        resolved_manifest_file,
        na.strings = c("", "NA", "N/A", "NULL", "."),
        check.names = FALSE
      ),
      traits = manifest[["trait"]],
      preparation_dir = preparation_dir,
      signature = signature,
      bcftools = provenance[["bcftools"]],
      bcftools_version = provenance[["bcftools_version"]],
      failed_alignment_file = normalizePath(failed_file, mustWork = TRUE),
      reused = TRUE
    ))
  }

  vcf_preparation <- prepare_manifest_vcf_files(
    params, manifest, preparation_dir
  )
  manifest <- vcf_preparation$manifest
  built <- build_manifest_fastasset_input(
    manifest, params, analysis_file, audit_file, allele_file, failed_file
  )
  data.table::fwrite(
    manifest, resolved_manifest_file,
    sep = "\t", quote = FALSE, na = "NA"
  )
  write_conversion_metadata(
    list(
      preparation_signature = signature,
      manifest = params$sumstats_manifest,
      number_traits = nrow(manifest),
      number_snps = built$n_snps,
      trait_order = paste(manifest[["trait"]], collapse = ";"),
      fastasset_neff_quantitative = "NEF/total analyzed N used directly",
      fastasset_neff_binary = "NCASE*NCONTROL/(NCASE+NCONTROL)",
      allele_policy = paste(
        "exact pair retained; swapped pair flips BETA; incompatible",
        "SNP-trait pair has BETA/SE/NEF set to NA and is logged"
      ),
      number_failed_allele_alignments = built$n_failed_alignments,
      failed_allele_alignment_file = built$failed_alignment_file,
      bcftools = vcf_preparation$bcftools,
      bcftools_version = vcf_preparation$bcftools_version,
      analysis_input = built$path
    ),
    provenance_file
  )
  writeLines(signature, completion_file, useBytes = TRUE)
  message("Constructed manifest-derived FastASSET input: ", built$path)

  list(
    path = built$path,
    manifest = manifest,
    traits = manifest[["trait"]],
    preparation_dir = preparation_dir,
    signature = signature,
    bcftools = vcf_preparation$bcftools,
    bcftools_version = vcf_preparation$bcftools_version,
    failed_alignment_file = built$failed_alignment_file,
    reused = FALSE
  )
}
