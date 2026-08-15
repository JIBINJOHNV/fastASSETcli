test_that("omitting ldsc-rdata selects validated automatic mode", {
  expect_error(
    fastASSETcli:::parse_fastasset_cli(c(
      "--sumstats-manifest", "manifest.tsv",
      "--output-dir", "output",
      "--run-name", "test"
    )),
    "automatic-LDSC mode requires"
  )
})

test_that("existing and automatic LDSC inputs cannot be mixed", {
  expect_error(
    fastASSETcli:::parse_fastasset_cli(c(
      "--sumstats-manifest", "manifest.tsv",
      "--ldsc-rdata", "ldsc.RData",
      "--hm3", "hm3.txt",
      "--ld-ref", "ld",
      "--output-dir", "output",
      "--run-name", "test"
    )),
    "cannot be combined"
  )
})

test_that("the removed fastasset-input option is rejected", {
  expect_error(
    fastASSETcli:::parse_fastasset_cli(c(
      "--fastasset-input", "wide.tsv.gz",
      "--sumstats-manifest", "manifest.tsv",
      "--output-dir", "output",
      "--run-name", "test"
    )),
    "Unknown option.*--fastasset-input"
  )
  expect_false(grepl(
    "--fastasset-input", fastASSETcli:::fastasset_cli_help(), fixed = TRUE
  ))
})

test_that("existing LDSC mode is fully specified by manifest and LDSC object", {
  directory <- tempfile("manifest_existing_ldsc_test_")
  dir.create(directory)
  summary_file <- file.path(directory, "trait.tsv")
  manifest <- file.path(directory, "manifest.tsv")
  ldsc <- file.path(directory, "ldsc.RData")
  writeLines("ID\tA1\tA2\tBETA\tSE\tNEF\nrs1\tA\tG\t0.1\t0.02\t100", summary_file)
  writeLines(
    "TRAIT\tFILE\ntrait_one\ttrait.tsv",
    manifest
  )
  writeLines("placeholder", ldsc)
  params <- fastASSETcli:::parse_fastasset_cli(c(
    "--sumstats-manifest", manifest,
    "--ldsc-rdata", ldsc,
    "--output-dir", directory,
    "--run-name", "test"
  ))
  expect_identical(params$ldsc_mode, "existing")
  expect_false("fastasset_input" %in% names(params))
})

test_that("GenomicSEM column mapping is parsed exactly", {
  mapping <- fastASSETcli:::parse_munge_column_map(
    "SNP=ID,A1=EA,A2=NEA,effect=BETA,SE=SE,N=NEF"
  )
  expect_identical(
    mapping,
    list(
      SNP = "ID", A1 = "EA", A2 = "NEA", effect = "BETA", SE = "SE",
      N = "NEF"
    )
  )
  expect_error(
    fastASSETcli:::parse_munge_column_map("N=NEF,N=SS"),
    "duplicate"
  )
})

test_that("manifest row order is canonical and explicit reordering remains validated", {
  directory <- tempfile("ldsc_manifest_test_")
  dir.create(directory)
  file_a <- file.path(directory, "a.tsv")
  file_b <- file.path(directory, "b.tsv")
  writeLines("SNP\tA1\tA2\tN\tZ\nrs1\tA\tG\t100\t1", file_a)
  writeLines("SNP\tA1\tA2\tN\tZ\nrs1\tA\tG\t100\t1", file_b)
  manifest <- file.path(directory, "manifest.tsv")
  writeLines(
    c(
      "NAME\tFILE\tSPREV\tPPREV",
      "trait_b\tb.tsv\tNA\tNA",
      "trait_a\ta.tsv\tNA\tNA"
    ),
    manifest
  )

  canonical <- fastASSETcli:::read_ldsc_manifest(manifest)
  expect_equal(as.character(canonical$trait), c("trait_b", "trait_a"))
  expect_identical(canonical$order, 1:2)

  resolved <- fastASSETcli:::read_ldsc_manifest(manifest, c("trait_a", "trait_b"))
  expect_equal(as.character(resolved$trait), c("trait_a", "trait_b"))
  expect_true(all(is.na(resolved$sample_prev)))
  expect_true(all(is.na(resolved$population_prev)))
  expect_identical(resolved$order, 1:2)
})

test_that("manifest tables construct wide input and align swapped effects", {
  directory <- tempfile("manifest_wide_test_")
  dir.create(directory)
  file_a <- file.path(directory, "a.tsv")
  file_b <- file.path(directory, "b.tsv")
  writeLines(
    c(
      "ID\tEA\tNEA\tBETA\tSE\tNEF",
      "rs1\tA\tG\t0.2\t0.05\t100",
      "rs2\tC\tT\t-0.1\t0.04\t90"
    ),
    file_a
  )
  writeLines(
    c(
      "ID\tEA\tNEA\tBETA\tSE\tNEF",
      "rs1\tG\tA\t0.3\t0.06\t200",
      "rs3\tA\tC\t0.4\t0.07\t180"
    ),
    file_b
  )
  manifest_file <- file.path(directory, "manifest.tsv")
  writeLines(
    c(
      "TRAIT\tFILE\tN\tSAMPLE\tSAMPLE_PREV\tPOPULATION_PREV",
      "trait_a\ta.tsv\tNA\tNA\tNA\tNA",
      "trait_b\tb.tsv\tNA\tNA\tNA\tNA"
    ),
    manifest_file
  )
  manifest <- fastASSETcli:::read_ldsc_manifest(manifest_file)
  output <- file.path(directory, "wide.tsv.gz")
  audit <- file.path(directory, "audit.tsv")
  alleles <- file.path(directory, "alleles.tsv.gz")
  built <- fastASSETcli:::build_manifest_fastasset_input(
    manifest,
    list(munge_column_map = list()),
    output,
    audit,
    alleles
  )
  wide <- data.table::fread(built$path)
  expect_identical(
    names(wide),
    c(
      "ID", "trait_a.Beta", "trait_a.SE", "trait_a.NEF",
      "trait_b.Beta", "trait_b.SE", "trait_b.NEF"
    )
  )
  expect_equal(wide[ID == "rs1", trait_b.Beta], -0.3)
  expect_true(is.na(wide[ID == "rs2", trait_b.Beta]))
  expect_equal(wide[ID == "rs3", trait_b.NEF], 180)
  expect_equal(data.table::fread(audit)[trait == "trait_b", allele_flips_to_reference], 1)
})

test_that("binary FastASSET Neff follows the original reference definition", {
  directory <- tempfile("binary_neff_test_")
  dir.create(directory)
  summary_file <- file.path(directory, "binary.tsv")
  writeLines(
    c(
      "ID\tA1\tA2\tBETA\tSE\tNC\tNCO\tN",
      "rs1\tA\tG\t0.2\t0.05\t100\t900\t1000"
    ),
    summary_file
  )
  row <- data.table::data.table(
    order = 1L, trait = "binary", source_file = summary_file,
    file = summary_file, source_format = "TABLE", N = NA_real_,
    sample_prev = 0.1, population_prev = 0.02
  )
  result <- fastASSETcli:::read_manifest_trait_for_fastasset(
    row, list(munge_column_map = list())
  )
  expect_equal(result$table$NEF, 90)
  expect_match(
    result$audit$fastasset_neff_source,
    "NC\\*NCO/\\(NC\\+NCO\\)"
  )
})

test_that("VCF manifests allow automatic binary sample prevalence", {
  directory <- tempfile("ldsc_vcf_manifest_test_")
  dir.create(directory)
  vcf <- file.path(directory, "binary.vcf")
  writeLines("##fileformat=VCFv4.2", vcf)
  manifest <- file.path(directory, "manifest.tsv")
  writeLines(
    c(
      "TRAIT\tFILE\tPOPULATION_PREV",
      "binary_trait\tbinary.vcf\t0.02"
    ),
    manifest
  )

  resolved <- fastASSETcli:::read_ldsc_manifest(
    manifest, "binary_trait"
  )
  expect_identical(resolved$source_format, "VCF")
  expect_true(is.na(resolved$sample_prev))
  expect_equal(resolved$population_prev, 0.02)
  expect_error(
    {
      writeLines(
        c(
          "TRAIT\tFILE\tN\tPOPULATION_PREV",
          "binary_trait\tbinary.vcf\t1000\t0.02"
        ),
        manifest
      )
      fastASSETcli:::read_ldsc_manifest(manifest, "binary_trait")
    },
    "Do not provide manifest N"
  )
})

test_that("VCF conversion applies allele, P, N and prevalence rules", {
  query <- data.table::data.table(
    CHR = c("1", "1", "1"),
    POS = c(100, 200, 300),
    SNP = c("rs1", "rs2", "rs3"),
    REF = c("A", "C", "G"),
    ALT = c("G", "T", "A"),
    ES = c(0.2, -0.1, 0.3),
    SE = c(0.05, 0.04, 0.06),
    LP = c(2, 1, 400),
    AF = c(0.2, 0.8, 0.4),
    NEF = c(360, 324, 360),
    NC = c(100, 90, 100),
    NCO = c(900, 810, 900),
    SS = c(1000, 900, 1000),
    SI = c(0.99, 0.95, 0.98)
  )
  binary <- fastASSETcli:::standardize_vcf_query_table(
    query,
    binary = TRUE,
    supplied_sample_prev = NA_real_,
    has_info = TRUE,
    source_file = "binary.vcf"
  )
  expect_equal(binary$sample_prev, 0.1)
  expect_equal(binary$table$N, c(1000, 900, 1000))
  expect_equal(binary$table$P[1L], 0.01)
  expect_gt(binary$table$P[3L], 0)
  expect_equal(binary$p_values_capped, 1)
  expect_identical(binary$table$A1, query$ALT)
  expect_identical(binary$table$A2, query$REF)

  quantitative <- fastASSETcli:::standardize_vcf_query_table(
    query[, setdiff(names(query), c("NC", "NCO")), with = FALSE],
    binary = FALSE,
    supplied_sample_prev = NA_real_,
    has_info = TRUE,
    source_file = "quantitative.vcf"
  )
  expect_true(is.na(quantitative$sample_prev))
  expect_equal(quantitative$table$N, query$NEF)
})

test_that("bcftools creates an audited LDSC table from GWAS-VCF", {
  bcftools <- unname(Sys.which("bcftools"))
  skip_if(!nzchar(bcftools), "bcftools is not installed")
  directory <- tempfile("bcftools_vcf_test_")
  dir.create(directory)
  vcf <- file.path(directory, "binary.vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##contig=<ID=1,length=1000000>",
      "##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"Effect\">",
      "##FORMAT=<ID=SE,Number=A,Type=Float,Description=\"SE\">",
      "##FORMAT=<ID=LP,Number=A,Type=Float,Description=\"-log10 P\">",
      "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"AF\">",
      "##FORMAT=<ID=SI,Number=A,Type=Float,Description=\"INFO\">",
      "##FORMAT=<ID=SS,Number=A,Type=Integer,Description=\"Total N\">",
      "##FORMAT=<ID=NC,Number=A,Type=Integer,Description=\"Cases\">",
      "##FORMAT=<ID=NCO,Number=A,Type=Integer,Description=\"Controls\">",
      "##FORMAT=<ID=NEF,Number=A,Type=Integer,Description=\"Effective N\">",
      paste0(
        "##SAMPLE=<ID=trait_one,StudyType=CaseControl,",
        "TotalCases=100,TotalControls=900>"
      ),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttrait_one",
      paste0(
        "1\t100\trs1\tA\tG\t.\tPASS\t.\t",
        "ES:SE:LP:AF:SI:SS:NC:NCO:NEF\t",
        "0.2:0.05:2:0.2:0.99:1000:100:900:360"
      )
    ),
    vcf
  )
  output <- file.path(directory, "converted.tsv.gz")
  metadata <- file.path(directory, "converted.metadata.tsv")
  result <- fastASSETcli:::prepare_one_vcf_summary(
    source_file = vcf,
    requested_sample = NA_character_,
    population_prev = 0.02,
    supplied_sample_prev = NA_real_,
    output_file = output,
    metadata_file = metadata,
    bcftools = bcftools,
    bcftools_version = "test"
  )
  converted <- data.table::fread(output)
  expect_equal(result$sample_prev, 0.1)
  expect_identical(result$sample, "trait_one")
  expect_equal(converted$P, 0.01)
  expect_equal(converted$N, 1000)
  expect_identical(converted$A1, "G")
  expect_identical(converted$A2, "A")
  expect_true(file.exists(metadata))

  manifest_row <- data.table::data.table(
    order = 1L, trait = "trait_one", source_file = vcf,
    file = output, source_format = "VCF", N = NA_real_,
    sample_prev = result$sample_prev, population_prev = 0.02
  )
  analysis <- fastASSETcli:::read_manifest_trait_for_fastasset(
    manifest_row, list(munge_column_map = list())
  )
  expect_equal(analysis$table$NEF, 90)
})

test_that("generated LDSC matrices receive exact two-axis names", {
  object <- list(I = diag(2), S = diag(2))
  named <- fastASSETcli:::add_ldsc_dimnames(object, c("trait_a", "trait_b"))
  expect_identical(rownames(named$I), c("trait_a", "trait_b"))
  expect_identical(colnames(named$I), c("trait_a", "trait_b"))
  expect_identical(rownames(named$S), c("trait_a", "trait_b"))
  expect_identical(colnames(named$S), c("trait_a", "trait_b"))
})
