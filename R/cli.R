fastasset_cli_help <- function() {
  paste0(
    "fastASSETcli: reference-aligned fastASSET for 200-250 traits\n\n",
    "Usage:\n",
    "  Existing LDSC object:\n",
    "  fastasset \\\n",
    "    --sumstats-manifest MANIFEST.tsv \\\n",
    "    --ldsc-rdata GenomicSEM_LDSCoutput.RData \\\n",
    "    --output-dir DIRECTORY \\\n",
    "    --run-name NAME [options]\n\n",
    "  Generate LDSC automatically:\n",
    "  fastasset \\\n",
    "    --sumstats-manifest MANIFEST.tsv \\\n",
    "    --hm3 w_hm3.snplist \\\n",
    "    --ld-ref LD_SCORE_DIRECTORY \\\n",
    "    --output-dir DIRECTORY \\\n",
    "    --run-name NAME [options]\n\n",
    "Always required:\n",
    "  --sumstats-manifest FILE     One TSV/CSV manifest row per trait. FILE\n",
    "                                may be a summary table or GWAS-VCF/BCF.\n",
    "                                The manifest row order is the canonical\n",
    "                                FastASSET and LDSC trait order.\n",
    "  --output-dir DIRECTORY       Output directory\n",
    "  --run-name NAME              Safe output prefix\n\n",
    "Input-preparation options (both LDSC modes):\n",
    "  --bcftools COMMAND           bcftools executable for VCF inputs\n",
    "                                [bcftools]\n",
    "  --vcf-cores INT              Parallel VCF conversions [--ncores]\n",
    "  --munge-column-map MAP       Tabular-column mapping used to construct\n",
    "                                FastASSET input and by GenomicSEM, e.g.\n",
    "                                SNP=ID,A1=EA,A2=NEA,effect=BETA,SE=SE,N=NEF\n\n",
    "Existing-LDSC mode:\n",
    "  --ldsc-rdata FILE            GenomicSEM .RData/.rda/.rds file\n",
    "  --ldsc-object-name NAME      R object containing $I [LDSCoutput]\n",
    "                                Trait names are recovered from $S if needed,\n",
    "                                then reordered to manifest trait order.\n\n",
    "Automatic-LDSC mode (used when --ldsc-rdata is omitted):\n",
    "  --hm3 FILE                   Uncompressed HapMap3 SNP reference\n",
    "  --ld-ref DIRECTORY           GenomicSEM LD-score reference directory\n",
    "  --wld-ref DIRECTORY          Separate weight directory [--ld-ref]\n",
    "  --munge-info-filter NUMBER   GenomicSEM munge INFO threshold [0.9]\n",
    "  --munge-maf-filter NUMBER    GenomicSEM munge MAF threshold [0.01]\n",
    "  --munge-cores INT            Parallel munging workers [--ncores]\n",
    "  --ldsc-chr INT               Chromosomes 1..INT [22]\n",
    "  --ldsc-blocks INT            Requested jackknife blocks [200]\n",
    "  --ldsc-chisq-max AUTO|NUMBER Maximum SNP chi-square [AUTO]\n",
    "  Generated .RData, resolved manifest, logs and provenance are saved\n",
    "  in a separate run-specific LDSC directory under --output-dir.\n\n",
    "VCF rules:\n",
    "  ES is the ALT-allele effect (A1=ALT, A2=REF); P=10^(-LP).\n",
    "  With POPULATION_PREV, NC/NCO define LDSC N and SAMPLE_PREV;\n",
    "  FastASSET Neff is NC*NCO/(NC+NCO), matching the reference.\n",
    "  Without POPULATION_PREV, quantitative N is FORMAT/NEF directly.\n\n",
    "Scientific options:\n",
    "  --scr-pthr NUMBER            Pre-screening P threshold [0.05]\n",
    "  --max-traits-per-side INT    Screened-trait guard per direction [16]\n",
    "  --min-available-traits INT   Minimum valid traits per SNP [2]\n",
    "  --cor-thr NUMBER             Correlation-block threshold [0.2]\n",
    "  --meth-pval DLM|IS|B         ASSET P-value method [DLM]\n",
    "  --include-meta TRUE|FALSE    Write screened fixed-effect Meta [TRUE]\n",
    "  --eigen-tolerance NUMBER     Positive-definite tolerance [1e-8]\n\n",
    "Compute and failure options:\n",
    "  --ncores INT                 Parallel SNP/chunk workers [90]\n",
    "  --chunk-size INT             SNP rows per chunk [1000]\n",
    "  --fail-on-critical TRUE|FALSE  Exit nonzero for CRITICAL QC [TRUE]\n",
    "  --fail-on-high TRUE|FALSE      Exit nonzero for HIGH QC [FALSE]\n",
    "  -h, --help                   Show this help\n",
    "  --version                    Show the installed package version\n\n",
    "Examples:\n",
    "  # Reuse an existing LDSC object\n",
    "  fastasset --sumstats-manifest /data/manifest.tsv \\\n",
    "    --ldsc-rdata /data/all_traits_LDSCoutput.RData \\\n",
    "    --output-dir /results/fastasset --run-name traits250 \\\n",
    "    --ncores 90 --chunk-size 1000 --include-meta TRUE\n\n",
    "  # Or generate LDSC, save it separately, then run FastASSET\n",
    "  fastasset --sumstats-manifest /data/manifest.tsv \\\n",
    "    --hm3 /refs/w_hm3.snplist --ld-ref /refs/eur_w_ld_chr \\\n",
    "    --output-dir /results/fastasset --run-name traits250 \\\n",
    "    --ncores 90 --munge-cores 90 --vcf-cores 90 --include-meta TRUE\n"
  )
}

parse_raw_cli <- function(arguments) {
  if (any(arguments %in% c("-h", "--help"))) {
    cat(fastasset_cli_help())
    quit(status = 0L)
  }

  values <- list()
  index <- 1L
  while (index <= length(arguments)) {
    token <- arguments[index]
    if (!startsWith(token, "--")) {
      stop("Unexpected positional argument: ", token, call. = FALSE)
    }

    if (grepl("=", token, fixed = TRUE)) {
      option <- sub("=.*$", "", token)
      value <- sub("^[^=]*=", "", token)
    } else {
      option <- token
      if (index == length(arguments) || startsWith(arguments[index + 1L], "--")) {
        stop("Option requires a value: ", option, call. = FALSE)
      }
      index <- index + 1L
      value <- arguments[index]
    }

    key <- sub("^--", "", option)
    if (!is.null(values[[key]])) {
      stop("Option supplied more than once: ", option, call. = FALSE)
    }
    values[[key]] <- value
    index <- index + 1L
  }
  values
}

parse_boolean <- function(value, option) {
  normalized <- tolower(value)
  if (normalized %in% c("true", "t", "1", "yes")) return(TRUE)
  if (normalized %in% c("false", "f", "0", "no")) return(FALSE)
  stop(option, " must be TRUE or FALSE.", call. = FALSE)
}

parse_integer <- function(value, option) {
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || !is.finite(parsed) || parsed != as.integer(parsed)) {
    stop(option, " must be an integer.", call. = FALSE)
  }
  as.integer(parsed)
}

parse_number <- function(value, option) {
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || !is.finite(parsed)) {
    stop(option, " must be a finite number.", call. = FALSE)
  }
  parsed
}

parse_auto_number <- function(value, option) {
  if (toupper(value) == "AUTO") return(NA_real_)
  parse_number(value, option)
}

parse_munge_column_map <- function(value) {
  if (is.null(value) || !nzchar(value) || toupper(value) == "AUTO") {
    return(list())
  }

  entries <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (any(!nzchar(entries))) {
    stop("--munge-column-map contains an empty entry.", call. = FALSE)
  }
  pairs <- strsplit(entries, "=", fixed = TRUE)
  if (any(lengths(pairs) != 2L)) {
    stop(
      "--munge-column-map must use KEY=COLUMN entries separated by commas.",
      call. = FALSE
    )
  }

  keys <- vapply(pairs, function(x) trimws(x[1L]), character(1L))
  values <- vapply(pairs, function(x) trimws(x[2L]), character(1L))
  keys_upper <- toupper(keys)
  allowed <- c(
    "SNP", "A1", "A2", "EFFECT", "SE", "INFO", "P", "N", "MAF", "Z",
    "NCASE", "NCONTROL"
  )
  invalid <- setdiff(keys_upper, allowed)
  if (length(invalid) > 0L) {
    stop(
      "Unsupported --munge-column-map key(s): ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(!nzchar(values))) {
    stop("--munge-column-map column names cannot be empty.", call. = FALSE)
  }
  if (anyDuplicated(keys_upper)) {
    stop("--munge-column-map contains duplicate keys.", call. = FALSE)
  }

  canonical_keys <- ifelse(keys_upper == "EFFECT", "effect", keys_upper)
  stats::setNames(as.list(values), canonical_keys)
}

parse_fastasset_cli <- function(arguments) {
  raw <- parse_raw_cli(arguments)
  allowed <- c(
    "sumstats-manifest", "ldsc-rdata", "ldsc-object-name", "output-dir",
    "hm3", "ld-ref", "wld-ref",
    "munge-info-filter", "munge-maf-filter", "munge-cores",
    "munge-column-map", "bcftools", "vcf-cores",
    "ldsc-chr", "ldsc-blocks", "ldsc-chisq-max",
    "run-name", "scr-pthr", "max-traits-per-side",
    "min-available-traits", "cor-thr", "meth-pval", "include-meta",
    "eigen-tolerance", "ncores", "chunk-size", "fail-on-critical",
    "fail-on-high"
  )
  unknown <- setdiff(names(raw), allowed)
  if (length(unknown) > 0L) {
    stop(
      "Unknown option(s): ", paste0("--", unknown, collapse = ", "),
      "\nUse --help for the accepted CLI.", call. = FALSE
    )
  }

  required <- c("sumstats-manifest", "output-dir", "run-name")
  supplied <- vapply(
    required,
    function(name) !is.null(raw[[name]]) && nzchar(raw[[name]]),
    logical(1L)
  )
  missing <- required[!supplied]
  if (length(missing) > 0L) {
    stop(
      "Required argument(s) not provided:\n  ",
      paste0("--", missing, collapse = "\n  "),
      "\n\nUse --help for a copy-ready example.", call. = FALSE
    )
  }

  has_existing_ldsc <- !is.null(raw[["ldsc-rdata"]]) &&
    nzchar(raw[["ldsc-rdata"]])
  automatic_required <- c("hm3", "ld-ref")
  automatic_supplied <- vapply(
    automatic_required,
    function(name) !is.null(raw[[name]]) && nzchar(raw[[name]]),
    logical(1L)
  )
  automatic_options <- c(
    automatic_required, "wld-ref", "munge-info-filter", "munge-maf-filter",
    "munge-cores", "ldsc-chr", "ldsc-blocks", "ldsc-chisq-max"
  )
  any_automatic_option <- any(names(raw) %in% automatic_options)

  if (has_existing_ldsc && any_automatic_option) {
    stop(
      "--ldsc-rdata cannot be combined with automatic-LDSC options. ",
      "Choose existing-LDSC mode or automatic-LDSC mode.",
      call. = FALSE
    )
  }
  if (!has_existing_ldsc && any(!automatic_supplied)) {
    missing_automatic <- automatic_required[!automatic_supplied]
    stop(
      "--ldsc-rdata was not supplied, so automatic-LDSC mode requires:\n  ",
      paste0("--", missing_automatic, collapse = "\n  "),
      "\n\nUse --help for a copy-ready example.",
      call. = FALSE
    )
  }

  get_value <- function(name, default) {
    if (is.null(raw[[name]])) default else raw[[name]]
  }

  parsed_ncores <- parse_integer(get_value("ncores", "90"), "--ncores")
  parsed_munge_cores <- parse_integer(
    get_value("munge-cores", as.character(parsed_ncores)), "--munge-cores"
  )
  params <- list(
    ldsc_rdata = raw[["ldsc-rdata"]],
    ldsc_mode = if (has_existing_ldsc) "existing" else "generate",
    ldsc_object_name = get_value("ldsc-object-name", "LDSCoutput"),
    sumstats_manifest = raw[["sumstats-manifest"]],
    hm3 = raw[["hm3"]],
    ld_ref = raw[["ld-ref"]],
    wld_ref = raw[["wld-ref"]],
    munge_info_filter = parse_number(
      get_value("munge-info-filter", "0.9"), "--munge-info-filter"
    ),
    munge_maf_filter = parse_number(
      get_value("munge-maf-filter", "0.01"), "--munge-maf-filter"
    ),
    munge_cores = parsed_munge_cores,
    munge_column_map = parse_munge_column_map(
      get_value("munge-column-map", "AUTO")
    ),
    bcftools = get_value("bcftools", "bcftools"),
    vcf_cores = parse_integer(
      get_value("vcf-cores", as.character(parsed_munge_cores)), "--vcf-cores"
    ),
    ldsc_chr = parse_integer(get_value("ldsc-chr", "22"), "--ldsc-chr"),
    ldsc_blocks = parse_integer(
      get_value("ldsc-blocks", "200"), "--ldsc-blocks"
    ),
    ldsc_chisq_max = parse_auto_number(
      get_value("ldsc-chisq-max", "AUTO"), "--ldsc-chisq-max"
    ),
    output_dir = raw[["output-dir"]],
    run_name = raw[["run-name"]],
    scr_pthr = parse_number(get_value("scr-pthr", "0.05"), "--scr-pthr"),
    max_numtraits_per_side = parse_integer(
      get_value("max-traits-per-side", "16"), "--max-traits-per-side"
    ),
    min_available_traits = parse_integer(
      get_value("min-available-traits", "2"), "--min-available-traits"
    ),
    cor_thr = parse_number(get_value("cor-thr", "0.2"), "--cor-thr"),
    meth_pval = toupper(get_value("meth-pval", "DLM")),
    include_meta = parse_boolean(
      get_value("include-meta", "TRUE"), "--include-meta"
    ),
    eigen_tolerance = parse_number(
      get_value("eigen-tolerance", "1e-8"), "--eigen-tolerance"
    ),
    ncores = parsed_ncores,
    chunk_size = parse_integer(get_value("chunk-size", "1000"), "--chunk-size"),
    fail_on_critical = parse_boolean(
      get_value("fail-on-critical", "TRUE"), "--fail-on-critical"
    ),
    fail_on_high = parse_boolean(
      get_value("fail-on-high", "FALSE"), "--fail-on-high"
    )
  )

  validate_cli_params(params)
}

validate_cli_params <- function(params) {
  missing_files <- c()
  if (!file.exists(params$sumstats_manifest)) {
    missing_files <- c(
      missing_files,
      paste0("--sumstats-manifest: ", params$sumstats_manifest)
    )
  }
  if (params$ldsc_mode == "existing" && !file.exists(params$ldsc_rdata)) {
    missing_files <- c(
      missing_files,
      paste0("--ldsc-rdata: ", params$ldsc_rdata)
    )
  }
  if (params$ldsc_mode == "generate" && !file.exists(params$hm3)) {
    missing_files <- c(missing_files, paste0("--hm3: ", params$hm3))
  }
  if (params$ldsc_mode == "generate" && !dir.exists(params$ld_ref)) {
    missing_files <- c(missing_files, paste0("--ld-ref: ", params$ld_ref))
  }
  if (params$ldsc_mode == "generate" && !is.null(params$wld_ref) &&
      !dir.exists(params$wld_ref)) {
    missing_files <- c(missing_files, paste0("--wld-ref: ", params$wld_ref))
  }
  if (length(missing_files) > 0L) {
    stop(
      "Input file(s) or directory/directories do not exist:\n  ",
      paste(missing_files, collapse = "\n  "),
      call. = FALSE
    )
  }

  unreadable_files <- c()
  if (file.access(params$sumstats_manifest, 4L) != 0L) {
    unreadable_files <- c(
      unreadable_files,
      paste0("--sumstats-manifest: ", params$sumstats_manifest)
    )
  }
  if (params$ldsc_mode == "existing" &&
      file.access(params$ldsc_rdata, 4L) != 0L) {
    unreadable_files <- c(
      unreadable_files,
      paste0("--ldsc-rdata: ", params$ldsc_rdata)
    )
  }
  if (params$ldsc_mode == "generate") {
    generated_inputs <- c(
      "--hm3" = params$hm3,
      "--ld-ref" = params$ld_ref,
      "--wld-ref" = if (is.null(params$wld_ref)) params$ld_ref else params$wld_ref
    )
    generated_unreadable <- generated_inputs[
      file.access(unname(generated_inputs), 4L) != 0L
    ]
    if (length(generated_unreadable) > 0L) {
      unreadable_files <- c(
        unreadable_files,
        paste0(names(generated_unreadable), ": ", generated_unreadable)
      )
    }
  }
  if (length(unreadable_files) > 0L) {
    stop(
      "Input file(s) or directory/directories are not readable:\n  ",
      paste(unreadable_files, collapse = "\n  "),
      call. = FALSE
    )
  }

  params$sumstats_manifest <- normalizePath(
    params$sumstats_manifest, mustWork = TRUE
  )
  if (params$ldsc_mode == "existing") {
    params$ldsc_rdata <- normalizePath(params$ldsc_rdata, mustWork = TRUE)
  } else {
    params$hm3 <- normalizePath(params$hm3, mustWork = TRUE)
    params$ld_ref <- normalizePath(params$ld_ref, mustWork = TRUE)
    params$wld_ref <- normalizePath(
      if (is.null(params$wld_ref)) params$ld_ref else params$wld_ref,
      mustWork = TRUE
    )
  }
  if (!dir.exists(params$output_dir)) {
    dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  params$output_dir <- normalizePath(params$output_dir, mustWork = TRUE)
  if (file.access(params$output_dir, 2L) != 0L) {
    stop("--output-dir is not writable: ", params$output_dir, call. = FALSE)
  }

  if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", params$run_name)) {
    stop("--run-name contains unsupported characters.", call. = FALSE)
  }
  if (!nzchar(params$ldsc_object_name)) {
    stop("--ldsc-object-name cannot be empty.", call. = FALSE)
  }
  if (params$scr_pthr <= 0 || params$scr_pthr >= 1) {
    stop("--scr-pthr must be strictly between 0 and 1.", call. = FALSE)
  }
  if (params$cor_thr <= 0 || params$cor_thr >= 1) {
    stop("--cor-thr must be strictly between 0 and 1.", call. = FALSE)
  }
  if (params$eigen_tolerance <= 0) {
    stop("--eigen-tolerance must be greater than zero.", call. = FALSE)
  }
  if (params$munge_info_filter < 0 || params$munge_info_filter > 1) {
    stop("--munge-info-filter must be between 0 and 1.", call. = FALSE)
  }
  if (params$munge_maf_filter < 0 || params$munge_maf_filter > 0.5) {
    stop("--munge-maf-filter must be between 0 and 0.5.", call. = FALSE)
  }
  if (params$munge_cores < 1L) {
    stop("--munge-cores must be at least 1.", call. = FALSE)
  }
  if (!nzchar(trimws(params$bcftools))) {
    stop("--bcftools cannot be empty.", call. = FALSE)
  }
  if (params$vcf_cores < 1L) {
    stop("--vcf-cores must be at least 1.", call. = FALSE)
  }
  if (params$ldsc_chr < 1L || params$ldsc_chr > 100L) {
    stop("--ldsc-chr must be between 1 and 100.", call. = FALSE)
  }
  if (params$ldsc_blocks < 2L) {
    stop("--ldsc-blocks must be at least 2.", call. = FALSE)
  }
  if (!is.na(params$ldsc_chisq_max) && params$ldsc_chisq_max <= 0) {
    stop("--ldsc-chisq-max must be AUTO or greater than zero.", call. = FALSE)
  }
  if (params$max_numtraits_per_side < 1L ||
      params$max_numtraits_per_side > 30L) {
    stop("--max-traits-per-side must be between 1 and 30.", call. = FALSE)
  }
  if (params$max_numtraits_per_side > 16L) {
    warning("A per-direction maximum above 16 is a sensitivity setting.")
  }
  if (params$min_available_traits < 1L) {
    stop("--min-available-traits must be at least 1.", call. = FALSE)
  }
  if (params$ncores < 1L) stop("--ncores must be at least 1.", call. = FALSE)
  if (params$chunk_size < 1L) {
    stop("--chunk-size must be at least 1.", call. = FALSE)
  }
  if (!params$meth_pval %in% c("DLM", "IS", "B")) {
    stop("--meth-pval must be DLM, IS or B.", call. = FALSE)
  }
  if (params$meth_pval == "IS" && params$max_numtraits_per_side > 10L) {
    stop("IS is documented as feasible only through 10 traits per direction.", call. = FALSE)
  }

  params
}
