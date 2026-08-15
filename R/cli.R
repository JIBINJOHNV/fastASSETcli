fastasset_cli_help <- function() {
  paste0(
    "fastASSETcli: reference-aligned fastASSET for 200-250 quantitative traits\n\n",
    "Usage:\n",
    "  fastasset \\\n",
    "    --fastasset-input FILE.tsv.gz \\\n",
    "    --ldsc-rdata GenomicSEM_LDSCoutput.RData \\\n",
    "    --output-dir DIRECTORY \\\n",
    "    --run-name NAME [options]\n\n",
    "Required:\n",
    "  --fastasset-input FILE       ID + <trait>.Beta/.SE/.NEF wide table\n",
    "  --ldsc-rdata FILE            Separate GenomicSEM .RData/.rda/.rds file\n",
    "  --output-dir DIRECTORY       Output directory\n",
    "  --run-name NAME              Safe output prefix\n\n",
    "LDSC extraction:\n",
    "  --ldsc-object-name NAME      R object containing $I [LDSCoutput]\n",
    "                                Trait names are recovered from $S if needed,\n",
    "                                then reordered to .Beta input order.\n\n",
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
    "Example:\n",
    "  fastasset --fastasset-input /data/asset.tsv.gz \\\n",
    "    --ldsc-rdata /data/all_traits_LDSCoutput.RData \\\n",
    "    --output-dir /results/fastasset --run-name traits250 \\\n",
    "    --ncores 90 --chunk-size 1000 --include-meta TRUE\n"
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

parse_fastasset_cli <- function(arguments) {
  raw <- parse_raw_cli(arguments)
  allowed <- c(
    "fastasset-input", "ldsc-rdata", "ldsc-object-name", "output-dir",
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

  required <- c("fastasset-input", "ldsc-rdata", "output-dir", "run-name")
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

  get_value <- function(name, default) {
    if (is.null(raw[[name]])) default else raw[[name]]
  }

  params <- list(
    fastasset_input = raw[["fastasset-input"]],
    ldsc_rdata = raw[["ldsc-rdata"]],
    ldsc_object_name = get_value("ldsc-object-name", "LDSCoutput"),
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
    ncores = parse_integer(get_value("ncores", "90"), "--ncores"),
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
  if (!file.exists(params$fastasset_input)) {
    missing_files <- c(missing_files, paste0("--fastasset-input: ", params$fastasset_input))
  }
  if (!file.exists(params$ldsc_rdata)) {
    missing_files <- c(missing_files, paste0("--ldsc-rdata: ", params$ldsc_rdata))
  }
  if (length(missing_files) > 0L) {
    stop("Input file(s) do not exist:\n  ", paste(missing_files, collapse = "\n  "),
         call. = FALSE)
  }

  unreadable_files <- c()
  if (file.access(params$fastasset_input, 4L) != 0L) {
    unreadable_files <- c(
      unreadable_files,
      paste0("--fastasset-input: ", params$fastasset_input)
    )
  }
  if (file.access(params$ldsc_rdata, 4L) != 0L) {
    unreadable_files <- c(
      unreadable_files,
      paste0("--ldsc-rdata: ", params$ldsc_rdata)
    )
  }
  if (length(unreadable_files) > 0L) {
    stop(
      "Input file(s) are not readable:\n  ",
      paste(unreadable_files, collapse = "\n  "),
      call. = FALSE
    )
  }

  params$fastasset_input <- normalizePath(params$fastasset_input, mustWork = TRUE)
  params$ldsc_rdata <- normalizePath(params$ldsc_rdata, mustWork = TRUE)
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
  if (params$chunk_size < 1L) stop("--chunk-size must be at least 1.", call. = FALSE)
  if (!params$meth_pval %in% c("DLM", "IS", "B")) {
    stop("--meth-pval must be DLM, IS or B.", call. = FALSE)
  }
  if (params$meth_pval == "IS" && params$max_numtraits_per_side > 10L) {
    stop("IS is documented as feasible only through 10 traits per direction.", call. = FALSE)
  }

  params
}
