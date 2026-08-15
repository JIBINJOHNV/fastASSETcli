test_that("omitting ldsc-rdata selects validated automatic mode", {
  expect_error(
    fastASSETcli:::parse_fastasset_cli(c(
      "--fastasset-input", "input.tsv.gz",
      "--output-dir", "output",
      "--run-name", "test"
    )),
    "automatic-LDSC mode requires"
  )
})

test_that("existing and automatic LDSC inputs cannot be mixed", {
  expect_error(
    fastASSETcli:::parse_fastasset_cli(c(
      "--fastasset-input", "input.tsv.gz",
      "--ldsc-rdata", "ldsc.RData",
      "--sumstats-manifest", "manifest.tsv",
      "--hm3", "hm3.txt",
      "--ld-ref", "ld",
      "--output-dir", "output",
      "--run-name", "test"
    )),
    "cannot be combined"
  )
})

test_that("GenomicSEM column mapping is parsed exactly", {
  mapping <- fastASSETcli:::parse_munge_column_map(
    "SNP=ID,A1=EA,A2=NEA,effect=BETA,N=NEF"
  )
  expect_identical(
    mapping,
    list(SNP = "ID", A1 = "EA", A2 = "NEA", effect = "BETA", N = "NEF")
  )
  expect_error(
    fastASSETcli:::parse_munge_column_map("N=NEF,N=SS"),
    "duplicate"
  )
})

test_that("manifest traits are validated and reordered", {
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

  resolved <- fastASSETcli:::read_ldsc_manifest(
    manifest, c("trait_a", "trait_b")
  )
  expect_equal(as.character(resolved$trait), c("trait_a", "trait_b"))
  expect_true(all(is.na(resolved$sample_prev)))
  expect_true(all(is.na(resolved$population_prev)))
  expect_identical(resolved$order, 1:2)
})

test_that("generated LDSC matrices receive exact two-axis names", {
  object <- list(I = diag(2), S = diag(2))
  named <- fastASSETcli:::add_ldsc_dimnames(object, c("trait_a", "trait_b"))
  expect_identical(rownames(named$I), c("trait_a", "trait_b"))
  expect_identical(colnames(named$I), c("trait_a", "trait_b"))
  expect_identical(rownames(named$S), c("trait_a", "trait_b"))
  expect_identical(colnames(named$S), c("trait_a", "trait_b"))
})
