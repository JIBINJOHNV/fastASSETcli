check_required_packages <- function() {
  packages <- c("data.table", "ASSET")
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1L), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop(
      "Install required R package(s): ",
      paste(missing_packages, collapse = ", ")
    )
  }
}

valid_trait_label_vector <- function(labels, expected_length) {
  !is.null(labels) && length(labels) == expected_length &&
    !anyNA(labels) && all(nzchar(as.character(labels))) &&
    !anyDuplicated(as.character(labels))
}

extract_genomicsem_ldsc <- function(path, object_name, analysis_traits) {
  analysis_traits <- as.character(analysis_traits)
  if (length(analysis_traits) == 0L || anyNA(analysis_traits) ||
      any(!nzchar(analysis_traits)) || anyDuplicated(analysis_traits)) {
    stop("The fastASSET analysis trait names must be non-empty and unique.")
  }

  loaded_object_name <- object_name
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    object <- readRDS(path)
    loaded_object_name <- "RDS_object"
  } else {
    environment <- new.env(parent = emptyenv())
    loaded_names <- load(path, envir = environment)
    if (object_name %in% loaded_names) {
      object <- get(object_name, envir = environment, inherits = FALSE)
    } else if (length(loaded_names) == 1L) {
      loaded_object_name <- loaded_names[1L]
      object <- get(loaded_object_name, envir = environment, inherits = FALSE)
      warning(
        "Object '", object_name, "' was not present; using the only object: ",
        loaded_object_name
      )
    } else {
      stop(
        "Object '", object_name, "' was not found in --ldsc-rdata. Objects: ",
        paste(loaded_names, collapse = ", "),
        ". Supply the correct --ldsc-object-name."
      )
    }
  }

  if (is.matrix(object)) {
    intercept <- object
    row_labels <- rownames(intercept)
    column_labels <- colnames(intercept)
    if (valid_trait_label_vector(row_labels, nrow(intercept)) &&
        valid_trait_label_vector(column_labels, ncol(intercept)) &&
        !identical(as.character(row_labels), as.character(column_labels))) {
      stop(
        "The direct LDSC intercept matrix has different row and column names; ",
        "their axis mapping is ambiguous. Supply matching dimnames."
      )
    }
    candidate_sources <- list(
      I_matching_dimnames = if (
        identical(as.character(row_labels), as.character(column_labels))
      ) column_labels else NULL,
      I_colnames_only = if (is.null(row_labels)) column_labels else NULL,
      I_rownames_only = if (is.null(column_labels)) row_labels else NULL
    )
  } else {
    if (!is.list(object) || is.null(object$I)) {
      stop(
        "The selected GenomicSEM object must be a list containing matrix I, ",
        "or it must itself be an intercept matrix."
      )
    }
    intercept <- as.matrix(object$I)
    # GenomicSEM::ldsc() constructs I and S in the same trait order, but
    # normally places trait.names only on S.  Therefore S column names are the
    # canonical label source for an unmodified GenomicSEM LDSC object.
    candidate_sources <- list(
      S_colnames = if (!is.null(object$S)) colnames(object$S) else NULL,
      S_rownames = if (!is.null(object$S)) rownames(object$S) else NULL,
      S_Stand_colnames = if (!is.null(object$S_Stand)) {
        colnames(object$S_Stand)
      } else NULL,
      S_Stand_rownames = if (!is.null(object$S_Stand)) {
        rownames(object$S_Stand)
      } else NULL,
      I_matching_dimnames = if (
        identical(rownames(intercept), colnames(intercept))
      ) colnames(intercept) else NULL,
      I_colnames = colnames(intercept),
      I_rownames = rownames(intercept)
    )
  }

  storage.mode(intercept) <- "double"
  if (nrow(intercept) != ncol(intercept)) {
    stop("GenomicSEM I must be a square matrix.")
  }
  matrix_size <- nrow(intercept)
  valid_candidates <- candidate_sources[vapply(
    candidate_sources,
    valid_trait_label_vector,
    logical(1L),
    expected_length = matrix_size
  )]
  if (length(valid_candidates) == 0L) {
    stop(
      "Could not recover unique trait names for GenomicSEM I. Checked I, S ",
      "and S_Stand dimnames."
    )
  }

  # Candidate sources are ordered by authority.  For a native GenomicSEM
  # object, S_colnames is used even if another object's fallback names would
  # happen to overlap more manifest labels: I and S share one native order.
  name_source <- names(valid_candidates)[1L]
  ldsc_traits <- as.character(valid_candidates[[1L]])
  rownames(intercept) <- colnames(intercept) <- ldsc_traits

  missing_from_ldsc <- setdiff(analysis_traits, ldsc_traits)
  missing_from_manifest <- setdiff(ldsc_traits, analysis_traits)
  common_traits <- analysis_traits[analysis_traits %in% ldsc_traits]

  if (length(missing_from_ldsc) > 0L) {
    warning(
      "Manifest trait(s) missing from LDSC and excluded: ",
      paste(missing_from_ldsc, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(missing_from_manifest) > 0L) {
    warning(
      "LDSC trait(s) missing from the manifest and excluded: ",
      paste(missing_from_manifest, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(common_traits) < 2L) {
    stop(
      "Fewer than two common traits remain between the manifest and LDSC ",
      "object."
    )
  }

  manifest_audit <- data.table::data.table(
    trait = analysis_traits,
    manifest_order = seq_along(analysis_traits),
    ldsc_order = match(analysis_traits, ldsc_traits),
    status = ifelse(
      analysis_traits %in% ldsc_traits, "COMMON", "MISSING_FROM_LDSC"
    ),
    analysis_order = match(analysis_traits, common_traits),
    action = ifelse(
      analysis_traits %in% ldsc_traits,
      "RETAINED_IN_MANIFEST_ORDER",
      "MANIFEST_TRAIT_EXCLUDED"
    )
  )
  ldsc_only_audit <- data.table::data.table(
    trait = missing_from_manifest,
    manifest_order = rep(NA_integer_, length(missing_from_manifest)),
    ldsc_order = match(missing_from_manifest, ldsc_traits),
    status = rep("MISSING_FROM_MANIFEST", length(missing_from_manifest)),
    analysis_order = rep(NA_integer_, length(missing_from_manifest)),
    action = rep("LDSC_TRAIT_EXCLUDED", length(missing_from_manifest))
  )
  trait_match_audit <- data.table::rbindlist(
    list(manifest_audit, ldsc_only_audit),
    use.names = TRUE,
    fill = TRUE
  )

  audit <- data.table::data.table(
    ldsc_original_order = seq_along(ldsc_traits),
    trait = ldsc_traits,
    used_in_fastasset = ldsc_traits %in% common_traits,
    fastasset_order = match(ldsc_traits, common_traits)
  )
  audit[, trait_name_source := name_source]
  audit[, ldsc_object_name := loaded_object_name]

  reordered <- intercept[common_traits, common_traits, drop = FALSE]
  stopifnot(identical(rownames(reordered), common_traits))
  stopifnot(identical(colnames(reordered), common_traits))

  message(
    "Extracted GenomicSEM I from object '", loaded_object_name,
    "'; trait names from ", name_source, ". Reordered ",
    length(common_traits), " common traits to manifest order."
  )

  list(
    intercept = reordered,
    order_audit = audit,
    trait_match_audit = trait_match_audit,
    common_traits = common_traits,
    missing_from_ldsc = missing_from_ldsc,
    missing_from_manifest = missing_from_manifest,
    object_name = loaded_object_name,
    trait_name_source = name_source,
    n_ldsc_traits = length(ldsc_traits),
    n_manifest_traits = length(analysis_traits),
    n_analysis_traits = length(common_traits)
  )
}

detect_fastasset_columns <- function(column_names) {
  beta_columns <- grep("\\.Beta$", column_names, value = TRUE)
  if (length(beta_columns) == 0L) {
    stop("No <trait>.Beta columns were found.")
  }

  traits <- sub("\\.Beta$", "", beta_columns)
  if (anyDuplicated(traits)) stop("Duplicate trait names were detected.")

  se_columns <- paste0(traits, ".SE")
  missing_se <- setdiff(se_columns, column_names)
  if (length(missing_se) > 0L) {
    stop("Missing SE columns: ", paste(missing_se, collapse = ", "))
  }

  nef_columns <- paste0(traits, ".NEF")
  legacy_n_columns <- paste0(traits, ".N")
  if (all(nef_columns %in% column_names)) {
    sample_size_columns <- nef_columns
    sample_size_suffix <- ".NEF"
  } else if (all(legacy_n_columns %in% column_names)) {
    sample_size_columns <- legacy_n_columns
    sample_size_suffix <- ".N"
    warning(
      "Using legacy <trait>.N columns as FastASSET NEF. Values are used ",
      "directly with no conversion. Rename them to <trait>.NEF for clarity."
    )
  } else {
    missing_nef <- setdiff(nef_columns, column_names)
    stop(
      "Missing FastASSET effective-sample-size columns. Expected ",
      "<trait>.NEF for every trait. Missing: ",
      paste(missing_nef, collapse = ", ")
    )
  }

  list(
    traits = traits,
    beta_columns = beta_columns,
    se_columns = se_columns,
    sample_size_columns = sample_size_columns,
    sample_size_suffix = sample_size_suffix
  )
}

empty_result_table <- function() {
  data.table::data.table(
    ID = character(), status = character(), severity = character(),
    screening_applied = logical(), analysis_scope = character(),
    p_two_sided = numeric(), p_positive = numeric(), p_negative = numeric(),
    beta_positive_adjusted = numeric(), se_positive_adjusted = numeric(),
    se_positive_meta = numeric(), beta_negative_adjusted = numeric(),
    se_negative_adjusted = numeric(), se_negative_meta = numeric(),
    n_selected_positive = integer(), n_selected_negative = integer(),
    selected_positive_traits = character(), selected_negative_traits = character(),
    n_screened_positive = integer(), n_screened_negative = integer(),
    screened_positive_traits = character(), screened_negative_traits = character(),
    candidate_subsets_total = numeric()
  )
}

empty_meta_table <- function() {
  data.table::data.table(
    ID = character(), status = character(), severity = character(),
    screening_applied = logical(), meta_scope = character(),
    n_analyzed_traits = integer(), analyzed_traits = character(),
    p_meta = numeric(), beta_meta_adjusted = numeric(),
    se_meta_adjusted = numeric(), n_screened_traits = integer(),
    n_screened_positive = integer(), n_screened_negative = integer(),
    screened_positive_traits = character(),
    screened_negative_traits = character()
  )
}

empty_qc_table <- function() {
  data.table::data.table(
    ID = character(), status = character(), severity = character(),
    message = character(), screening_applied = logical(),
    analysis_scope = character(), n_input_traits = integer(),
    n_valid_traits = integer(),
    n_invalid_traits = integer(), invalid_traits = character(),
    n_screened_traits = integer(), n_screened_positive = integer(),
    n_screened_negative = integer(), screened_positive_traits = character(),
    screened_negative_traits = character(), candidate_subsets_positive = numeric(),
    candidate_subsets_negative = numeric(), candidate_subsets_total = numeric(),
    chunk = integer(), runtime_seconds = numeric()
  )
}

make_direction_limit_table <- function(qc_table, max_traits_per_direction,
                                       scr_pthr) {
  table <- data.table::copy(
    qc_table[qc_table[["status"]] == "DIRECTION_LIMIT_EXCEEDED"]
  )
  table[["max_traits_per_direction"]] <- rep(
    as.integer(max_traits_per_direction), nrow(table)
  )
  table[["scr_pthr"]] <- rep(as.numeric(scr_pthr), nrow(table))
  table[["action"]] <- rep(
    "SNP_SKIPPED_ASSET_SEARCH_ANALYSIS_CONTINUED", nrow(table)
  )
  table
}

make_critical_qc <- function(snp, input_traits, valid_traits, message) {
  qc <- make_fastasset_qc(
    snp = snp,
    status = "ASSET_ERROR",
    severity = "CRITICAL",
    message = message,
    n_input_traits = length(input_traits),
    valid_traits = valid_traits,
    input_traits = input_traits
  )
  qc
}

run_fastasset_chunk <- function(chunk_id, chunk_input_file, output_dir,
                                columns, correlation, blocks, params) {
  result_file <- file.path(
    output_dir, paste0("fastasset_chunk_", chunk_id, "_results.tsv.gz")
  )
  qc_file <- file.path(
    output_dir, paste0("fastasset_chunk_", chunk_id, "_qc.tsv.gz")
  )
  meta_file <- file.path(
    output_dir, paste0("fastasset_chunk_", chunk_id, "_meta.tsv.gz")
  )
  direction_limit_file <- file.path(
    output_dir,
    paste0("fastasset_chunk_", chunk_id, "_direction_limit.tsv.gz")
  )
  done_file <- file.path(
    output_dir, paste0("fastasset_chunk_", chunk_id, ".done")
  )

  if (file.exists(done_file) && file.exists(result_file) && file.exists(qc_file) &&
      file.exists(meta_file) && file.exists(direction_limit_file)) {
    return(TRUE)
  }

  unlink(
    c(result_file, qc_file, meta_file, direction_limit_file, done_file),
    force = TRUE
  )
  chunk <- data.table::fread(chunk_input_file, check.names = FALSE)

  results <- vector("list", nrow(chunk))
  meta_results <- vector("list", nrow(chunk))
  qc_rows <- vector("list", nrow(chunk))
  result_index <- 0L
  meta_index <- 0L
  beta_columns <- columns$beta_columns
  se_columns <- columns$se_columns
  sample_size_columns <- columns$sample_size_columns

  for (row_index in seq_len(nrow(chunk))) {
    started <- proc.time()[["elapsed"]]
    snp <- as.character(chunk$ID[row_index])

    beta <- as.numeric(unlist(
      chunk[row_index, ..beta_columns],
      use.names = FALSE
    ))
    standard_error <- as.numeric(unlist(
      chunk[row_index, ..se_columns],
      use.names = FALSE
    ))
    # NEF was prepared according to the upstream trait-type-specific rule.
    nef <- as.numeric(unlist(
      chunk[row_index, ..sample_size_columns],
      use.names = FALSE
    ))
    valid_traits_for_qc <- columns$traits[
      is.finite(beta) & is.finite(standard_error) & is.finite(nef) &
        standard_error > 0 & nef > 0
    ]

    outcome <- tryCatch(
      fast_asset_2(
        snp = snp,
        traits.lab = columns$traits,
        beta.hat = beta,
        sigma.hat = standard_error,
        Neff = nef,
        cor = correlation,
        block = blocks,
        scr_pthr = params$scr_pthr,
        max_numtraits_per_side = params$max_numtraits_per_side,
        min_available_traits = params$min_available_traits,
        meth_pval = params$meth_pval,
        include_meta = params$include_meta
      ),
      error = function(error) {
        list(
          result = NULL,
          qc = make_critical_qc(
            snp, columns$traits, valid_traits_for_qc, error$message
          )
        )
      }
    )

    elapsed <- proc.time()[["elapsed"]] - started
    outcome$qc$chunk <- as.integer(chunk_id)
    outcome$qc$runtime_seconds <- as.numeric(elapsed)
    qc_rows[[row_index]] <- data.table::as.data.table(outcome$qc)

    extracted <- extract_fastasset_two_sided(outcome)
    if (!is.null(extracted)) {
      result_index <- result_index + 1L
      results[[result_index]] <- data.table::as.data.table(extracted)
    }

    extracted_meta <- extract_fastasset_meta(outcome)
    if (!is.null(extracted_meta)) {
      meta_index <- meta_index + 1L
      meta_results[[meta_index]] <- data.table::as.data.table(extracted_meta)
    }
  }

  qc_table <- if (length(qc_rows) > 0L) {
    data.table::rbindlist(qc_rows, fill = TRUE, use.names = TRUE)
  } else {
    empty_qc_table()
  }
  result_table <- if (result_index > 0L) {
    data.table::rbindlist(
      results[seq_len(result_index)], fill = TRUE, use.names = TRUE
    )
  } else {
    empty_result_table()
  }
  meta_table <- if (meta_index > 0L) {
    data.table::rbindlist(
      meta_results[seq_len(meta_index)], fill = TRUE, use.names = TRUE
    )
  } else {
    empty_meta_table()
  }
  direction_limit_table <- make_direction_limit_table(
    qc_table,
    max_traits_per_direction = params$max_numtraits_per_side,
    scr_pthr = params$scr_pthr
  )

  data.table::fwrite(qc_table, qc_file, sep = "\t", quote = FALSE, na = "NA")
  data.table::fwrite(
    result_table, result_file, sep = "\t", quote = FALSE, na = "NA"
  )
  data.table::fwrite(
    meta_table, meta_file, sep = "\t", quote = FALSE, na = "NA"
  )
  data.table::fwrite(
    direction_limit_table,
    direction_limit_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  writeLines(
    c(
      "status=complete",
      paste0("chunk=", chunk_id),
      paste0("input_rows=", nrow(chunk)),
      paste0("result_rows=", nrow(result_table)),
      paste0("meta_rows=", nrow(meta_table)),
      paste0("direction_limit_rows=", nrow(direction_limit_table))
    ),
    done_file
  )
  TRUE
}

summarize_qc_files <- function(qc_files) {
  counts <- lapply(qc_files[file.exists(qc_files)], function(file) {
    qc <- data.table::fread(file, select = c("status", "severity"))
    qc[, .N, by = .(status, severity)]
  })
  if (length(counts) == 0L) {
    return(data.table::data.table(
      status = character(), severity = character(), N = integer()
    ))
  }
  data.table::rbindlist(counts)[, .(N = sum(N)), by = .(status, severity)]
}

run_fastasset_stage <- function(analysis_input, correlation, blocks,
                                output_dir, columns, params) {
  if (.Platform$OS.type == "windows" && params$ncores > 1L) {
    stop("This parallel implementation uses forked workers and requires Linux/macOS.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  input_chunk_dir <- file.path(output_dir, "input_chunks")
  result_chunk_dir <- file.path(output_dir, "result_chunks")
  dir.create(input_chunk_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_chunk_dir, recursive = TRUE, showWarnings = FALSE)

  message("Reading manifest-derived FastASSET input: ", analysis_input)
  input <- data.table::fread(analysis_input, check.names = FALSE)
  if (!"ID" %in% names(input)) stop("Prepared analysis input is missing the ID column.")
  if (anyNA(input$ID) || any(!nzchar(as.character(input$ID)))) {
    stop("Every ID value must be non-missing and non-empty.")
  }
  if (anyDuplicated(input$ID)) {
    warning("Duplicate ID values exist; each row will still be analyzed independently.")
  }

  row_groups <- split(
    seq_len(nrow(input)),
    ceiling(seq_len(nrow(input)) / params$chunk_size)
  )
  n_chunks <- length(row_groups)
  if (n_chunks == 0L) stop("Prepared analysis input has no rows.")

  message(
    "Materializing ", n_chunks, " input chunks for ", nrow(input),
    " SNP rows and ", length(columns$traits), " traits."
  )
  chunk_input_files <- character(n_chunks)
  required_columns <- c(
    "ID", columns$beta_columns, columns$se_columns,
    columns$sample_size_columns
  )

  for (chunk_id in seq_len(n_chunks)) {
    chunk_input_files[chunk_id] <- file.path(
      input_chunk_dir, paste0("fastasset_input_chunk_", chunk_id, ".tsv.gz")
    )
    if (!file.exists(chunk_input_files[chunk_id])) {
      data.table::fwrite(
        input[row_groups[[chunk_id]], ..required_columns],
        chunk_input_files[chunk_id],
        sep = "\t",
        quote = FALSE,
        na = "NA"
      )
    }
  }
  rm(input, row_groups)
  invisible(gc())

  workers <- min(params$ncores, n_chunks)
  message(
    "Running ", workers, " parallel chunk workers. ",
    "Subset enumeration within each SNP remains sequential."
  )

  chunk_results <- parallel::mclapply(
    seq_len(n_chunks),
    function(chunk_id) {
      tryCatch(
        run_fastasset_chunk(
          chunk_id = chunk_id,
          chunk_input_file = chunk_input_files[chunk_id],
          output_dir = result_chunk_dir,
          columns = columns,
          correlation = correlation,
          blocks = blocks,
          params = params
        ),
        error = function(error) {
          message("Chunk ", chunk_id, " failed: ", error$message)
          FALSE
        }
      )
    },
    mc.cores = workers,
    mc.preschedule = FALSE,
    mc.allow.recursive = FALSE
  )

  failed_chunks <- which(!vapply(chunk_results, isTRUE, logical(1L)))
  if (length(failed_chunks) > 0L) {
    stop(
      "Chunk failure(s): ", paste(failed_chunks, collapse = ", "),
      ". Completed chunks are preserved and will be reused on rerun."
    )
  }

  result_chunk_files <- file.path(
    result_chunk_dir,
    paste0("fastasset_chunk_", seq_len(n_chunks), "_results.tsv.gz")
  )
  qc_chunk_files <- file.path(
    result_chunk_dir,
    paste0("fastasset_chunk_", seq_len(n_chunks), "_qc.tsv.gz")
  )
  meta_chunk_files <- file.path(
    result_chunk_dir,
    paste0("fastasset_chunk_", seq_len(n_chunks), "_meta.tsv.gz")
  )
  direction_limit_chunk_files <- file.path(
    result_chunk_dir,
    paste0(
      "fastasset_chunk_", seq_len(n_chunks), "_direction_limit.tsv.gz"
    )
  )

  final_result_file <- file.path(
    output_dir, paste0(params$run_name, "_fastasset_results.tsv.gz")
  )
  final_qc_file <- file.path(
    output_dir, paste0(params$run_name, "_fastasset_qc.tsv.gz")
  )
  final_meta_file <- file.path(
    output_dir, paste0(params$run_name, "_fastasset_meta.tsv.gz")
  )
  final_direction_limit_file <- file.path(
    output_dir,
    paste0(params$run_name, "_direction_limit_exceeded.tsv.gz")
  )

  combine_tsv_gz(result_chunk_files, final_result_file)
  combine_tsv_gz(qc_chunk_files, final_qc_file)
  combine_tsv_gz(meta_chunk_files, final_meta_file)
  combine_tsv_gz(
    direction_limit_chunk_files,
    final_direction_limit_file
  )

  summary <- summarize_qc_files(qc_chunk_files)
  summary_file <- file.path(
    output_dir, paste0(params$run_name, "_fastasset_summary.tsv")
  )
  data.table::fwrite(summary, summary_file, sep = "\t")

  critical_count <- summary[severity == "CRITICAL", sum(N)]
  high_count <- summary[severity == "HIGH", sum(N)]
  if (length(critical_count) == 0L || is.na(critical_count)) critical_count <- 0L
  if (length(high_count) == 0L || is.na(high_count)) high_count <- 0L

  if (high_count > 0L) {
    warning(
      high_count,
      " SNP row(s) exceeded the per-direction screening limit and were not tested. ",
      "The analysis continued; review ", final_direction_limit_file, "."
    )
  }

  outputs <- list(
    results = final_result_file,
    meta = final_meta_file,
    qc = final_qc_file,
    direction_limit_exceeded = final_direction_limit_file,
    summary = summary_file,
    output_dir = output_dir,
    critical_count = critical_count,
    high_count = high_count
  )

  if (params$fail_on_critical && critical_count > 0L) {
    stop(
      critical_count,
      " CRITICAL SNP-level error(s) were recorded. Outputs were written; ",
      "inspect ", final_qc_file
    )
  }
  if (params$fail_on_high && high_count > 0L) {
    stop(
      high_count,
      " HIGH-severity direction-limit exclusion(s) were recorded. Outputs were ",
      "written; inspect ", final_direction_limit_file
    )
  }

  outputs
}

run_asset_pipeline <- function(params) {
  check_required_packages()

  # In existing-LDSC mode, determine the common trait set before reading or
  # converting summary statistics.  Manifest order remains canonical, and
  # manifest rows (including VCF SAMPLE values) not represented in LDSC are
  # excluded from all downstream preparation.
  manifest_traits <- read_ldsc_manifest_trait_names(params$sumstats_manifest)
  trait_match_file <- file.path(
    params$output_dir,
    paste0(params$run_name, "_manifest_ldsc_trait_match.tsv")
  )
  ldsc_extraction <- NULL
  if (identical(params$ldsc_mode, "existing")) {
    ldsc_extraction <- extract_genomicsem_ldsc(
      path = params$ldsc_rdata,
      object_name = params$ldsc_object_name,
      analysis_traits = manifest_traits
    )
    data.table::fwrite(
      ldsc_extraction$trait_match_audit,
      trait_match_file,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )
    manifest <- read_ldsc_manifest(
      params$sumstats_manifest,
      retain_traits = ldsc_extraction$common_traits
    )
  } else {
    manifest <- read_ldsc_manifest(params$sumstats_manifest)
  }

  prepared_input <- prepare_manifest_analysis_input(params, manifest = manifest)
  params$analysis_input <- prepared_input$path
  params$input_preparation_dir <- prepared_input$preparation_dir
  params$input_preparation_signature <- prepared_input$signature
  params$variant_map_file <- prepared_input$variant_map_file
  params$failed_alignment_file <- prepared_input$failed_alignment_file

  input_header <- names(data.table::fread(
    params$analysis_input, nrows = 0L, check.names = FALSE
  ))
  columns <- detect_fastasset_columns(input_header)
  if (!identical(columns$traits, prepared_input$traits)) {
    stop(
      "Prepared FastASSET trait order does not exactly match manifest row order."
    )
  }
  if (length(columns$traits) < params$min_available_traits) {
    stop(
      "Input has only ", length(columns$traits), " traits but min_available_traits is ",
      params$min_available_traits, "."
    )
  }

  message(
    "Detected ", length(columns$traits), " manifest traits; sample-size suffix is ",
    columns$sample_size_suffix,
    ". Prepared NEF values follow the upstream fastASSET definition."
  )

  ldsc_source <- resolve_ldsc_source(
    params, columns$traits, prepared_input
  )
  params$ldsc_rdata <- ldsc_source$path
  params$ldsc_generation_dir <- ldsc_source$generation_dir
  params$ldsc_generation_signature <- ldsc_source$signature

  if (is.null(ldsc_extraction)) {
    ldsc_extraction <- extract_genomicsem_ldsc(
      path = params$ldsc_rdata,
      object_name = params$ldsc_object_name,
      analysis_traits = columns$traits
    )
    data.table::fwrite(
      ldsc_extraction$trait_match_audit,
      trait_match_file,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )
  }
  if (!identical(ldsc_extraction$common_traits, columns$traits)) {
    stop(
      "Internal trait-order mismatch after manifest/LDSC common-trait ",
      "selection."
    )
  }

  params$trait_match_file <- trait_match_file

  correlation <- normalize_ldsc_intercept(
    ldsc_extraction$intercept,
    traits = columns$traits,
    eigen_tolerance = params$eigen_tolerance
  )
  blocks <- create_blocks(correlation, cor_thr = params$cor_thr)

  signature <- make_run_signature(params, columns$traits, correlation)
  analysis_dir <- file.path(
    params$output_dir, paste0(params$run_name, "_", signature)
  )
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(
    try(capture_environment(analysis_dir), silent = TRUE),
    add = TRUE
  )

  normalized_correlation_file <- file.path(
    analysis_dir, "trait_correlation_normalized.tsv.gz"
  )
  reordered_intercept_file <- file.path(
    analysis_dir, "ldsc_intercept_I_reordered.tsv.gz"
  )
  ldsc_order_file <- file.path(
    analysis_dir, "ldsc_trait_order_audit.tsv"
  )
  analysis_trait_match_file <- file.path(
    analysis_dir, "manifest_ldsc_trait_match.tsv"
  )
  reordered_intercept_table <- data.table::as.data.table(
    ldsc_extraction$intercept,
    keep.rownames = "trait"
  )
  data.table::fwrite(
    reordered_intercept_table,
    reordered_intercept_file,
    sep = "\t",
    quote = FALSE
  )
  data.table::fwrite(
    ldsc_extraction$order_audit,
    ldsc_order_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  data.table::fwrite(
    ldsc_extraction$trait_match_audit,
    analysis_trait_match_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  correlation_table <- data.table::as.data.table(correlation, keep.rownames = "trait")
  data.table::fwrite(
    correlation_table,
    normalized_correlation_file,
    sep = "\t",
    quote = FALSE
  )
  saveRDS(blocks, file.path(analysis_dir, "trait_correlation_blocks.rds"))

  settings <- data.table::data.table(
    parameter = c(
      "run_signature", "number_input_traits", "sample_size_suffix",
      "scr_pthr", "max_numtraits_per_side", "min_available_traits",
      "cor_thr", "meth_pval", "include_meta", "ncores", "chunk_size",
      "minimum_correlation_eigenvalue", "ldsc_mode", "ldsc_rdata",
      "ldsc_object_name", "ldsc_generation_dir",
      "ldsc_generation_signature", "ldsc_trait_name_source",
      "number_manifest_traits", "number_ldsc_traits",
      "number_common_traits", "sumstats_manifest", "analysis_input",
      "input_preparation_dir", "input_preparation_signature",
      "variant_map_file", "failed_alignment_file",
      "manifest_ldsc_trait_match_file"
    ),
    value = as.character(c(
      signature, length(columns$traits), columns$sample_size_suffix,
      params$scr_pthr, params$max_numtraits_per_side,
      params$min_available_traits, params$cor_thr, params$meth_pval,
      params$include_meta, params$ncores, params$chunk_size,
      attr(correlation, "minimum_eigenvalue"), params$ldsc_mode,
      params$ldsc_rdata, ldsc_extraction$object_name,
      params$ldsc_generation_dir, params$ldsc_generation_signature,
      ldsc_extraction$trait_name_source, ldsc_extraction$n_manifest_traits,
      ldsc_extraction$n_ldsc_traits, ldsc_extraction$n_analysis_traits,
      params$sumstats_manifest, params$analysis_input,
      params$input_preparation_dir, params$input_preparation_signature,
      params$variant_map_file, params$failed_alignment_file,
      analysis_trait_match_file
    ))
  )
  data.table::fwrite(
    settings,
    file.path(analysis_dir, "analysis_settings.tsv"),
    sep = "\t"
  )

  outputs <- run_fastasset_stage(
    analysis_input = params$analysis_input,
    correlation = correlation,
    blocks = blocks,
    output_dir = analysis_dir,
    columns = columns,
    params = params
  )
  outputs$ldsc_rdata <- params$ldsc_rdata
  outputs$ldsc_generation_dir <- params$ldsc_generation_dir
  outputs$analysis_input <- params$analysis_input
  outputs$input_preparation_dir <- params$input_preparation_dir
  outputs$variant_map_file <- params$variant_map_file
  outputs$failed_alignment_file <- params$failed_alignment_file
  outputs$trait_match_file <- analysis_trait_match_file

  message("fastASSET analysis completed: ", analysis_dir)
  invisible(outputs)
}
