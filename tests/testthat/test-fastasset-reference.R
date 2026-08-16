identity_correlation <- function(traits) {
  correlation <- diag(length(traits))
  dimnames(correlation) <- list(traits, traits)
  correlation
}

test_that("one valid trait follows the published unscreened closed form", {
  skip_if_not_installed("ASSET")
  traits <- c("trait_a", "trait_b")
  outcome <- fastASSETcli:::fast_asset_2(
    snp = "rs_single_valid",
    traits.lab = traits,
    beta.hat = c(0.1, NA_real_),
    sigma.hat = c(0.05, NA_real_),
    Neff = c(100, NA_real_),
    cor = identity_correlation(traits),
    block = list(traits),
    min_available_traits = 1L,
    include_meta = TRUE
  )

  expected_p <- 2 * stats::pnorm(-2)
  extracted <- fastASSETcli:::extract_fastasset_two_sided(outcome)
  extracted_meta <- fastASSETcli:::extract_fastasset_meta(outcome)

  expect_identical(outcome$qc$status, "SINGLE_VALID_TRAIT_NO_SCREEN")
  expect_identical(outcome$qc$severity, "OK")
  expect_false(outcome$qc$screening_applied)
  expect_identical(
    outcome$qc$analysis_scope, "SINGLE_VALID_TRAIT_UNSCREENED"
  )
  expect_equal(extracted$p_two_sided, expected_p)
  expect_equal(extracted$p_positive, expected_p)
  expect_equal(extracted$p_negative, 1)
  expect_identical(extracted$selected_positive_traits, "trait_a")
  expect_identical(extracted$selected_negative_traits, "")
  expect_equal(extracted_meta$p_meta, expected_p)
  expect_equal(extracted_meta$beta_meta_adjusted, 0.2)
  expect_equal(extracted_meta$se_meta_adjusted, 0.1)
  expect_identical(extracted_meta$n_analyzed_traits, 1L)
  expect_identical(extracted_meta$analyzed_traits, "trait_a")
  expect_identical(
    extracted_meta$meta_scope, "SINGLE_VALID_TRAIT_UNSCREENED"
  )
})

test_that("one screened trait uses its adjusted closed-form z statistic", {
  skip_if_not_installed("ASSET")
  traits <- c("trait_negative", "trait_null")
  outcome <- fastASSETcli:::fast_asset_2(
    snp = "rs_single_screened",
    traits.lab = traits,
    beta.hat = c(-0.3, 0.01),
    sigma.hat = c(0.1, 0.1),
    Neff = c(100, 100),
    cor = identity_correlation(traits),
    block = list(),
    min_available_traits = 1L,
    include_meta = TRUE
  )

  extracted <- fastASSETcli:::extract_fastasset_two_sided(outcome)
  adjusted_beta <- outcome$result$beta.hat[["trait_negative"]]
  adjusted_se <- outcome$result$sigma.hat[["trait_negative"]]
  expected_p <- 2 * stats::pnorm(-abs(adjusted_beta / adjusted_se))

  expect_identical(outcome$qc$status, "SINGLE_TRAIT_AFTER_SCREEN")
  expect_true(outcome$qc$screening_applied)
  expect_identical(
    outcome$qc$analysis_scope,
    "SINGLE_SCREENED_SELECTION_ADJUSTED_TRAIT"
  )
  expect_equal(extracted$p_two_sided, expected_p)
  expect_equal(extracted$p_positive, 1)
  expect_equal(extracted$p_negative, expected_p)
  expect_identical(extracted$selected_positive_traits, "")
  expect_identical(extracted$selected_negative_traits, "trait_negative")
})

test_that("zero valid traits and a stricter user minimum remain auditable", {
  skip_if_not_installed("ASSET")
  traits <- c("trait_a", "trait_b")
  correlation <- identity_correlation(traits)

  none <- fastASSETcli:::fast_asset_2(
    snp = "rs_none",
    traits.lab = traits,
    beta.hat = c(NA_real_, Inf),
    sigma.hat = c(NA_real_, 0),
    Neff = c(NA_real_, -1),
    cor = correlation,
    block = list()
  )
  expect_identical(none$qc$status, "NO_VALID_TRAITS")
  expect_null(none$result)

  strict <- fastASSETcli:::fast_asset_2(
    snp = "rs_strict",
    traits.lab = traits,
    beta.hat = c(0.1, NA_real_),
    sigma.hat = c(0.05, NA_real_),
    Neff = c(100, NA_real_),
    cor = correlation,
    block = list(),
    min_available_traits = 2L
  )
  expect_identical(strict$qc$status, "INSUFFICIENT_VALID_TRAITS")
  expect_null(strict$result)
})

test_that("zero adjusted effect belongs to the h.traits2 non-negative side", {
  result <- fastASSETcli:::make_single_trait_fastasset_result(
    snp = "rs_zero", trait = "trait_zero", beta = 0, sigma = 0.1
  )
  outcome <- list(
    result = result,
    qc = fastASSETcli:::make_fastasset_qc(
      "rs_zero", "SINGLE_VALID_TRAIT_NO_SCREEN", "OK", "test",
      1L, "trait_zero", "trait_zero",
      screening_applied = FALSE,
      analysis_scope = "SINGLE_VALID_TRAIT_UNSCREENED"
    )
  )
  extracted <- fastASSETcli:::extract_fastasset_two_sided(outcome)

  expect_equal(extracted$p_two_sided, 1)
  expect_equal(extracted$p_positive, 1)
  expect_equal(extracted$p_negative, 1)
  expect_identical(extracted$selected_positive_traits, "trait_zero")
  expect_identical(extracted$selected_negative_traits, "")
})

test_that("direction guard permits 16 and rejects 17 traits per side", {
  expect_false(fastASSETcli:::direction_limit_exceeded(rep(1, 16), 16L))
  expect_true(fastASSETcli:::direction_limit_exceeded(rep(1, 17), 16L))
  expect_false(fastASSETcli:::direction_limit_exceeded(rep(-1, 16), 16L))
  expect_true(fastASSETcli:::direction_limit_exceeded(rep(-1, 17), 16L))
})

test_that("independent screening protects extreme z statistics from underflow", {
  screened <- fastASSETcli:::scr_transform(
    betahat = c(extreme = 40),
    SE = c(extreme = 1),
    SNP = "rs_extreme",
    traits_i = "extreme",
    cor = identity_correlation("extreme"),
    block = list(),
    Neff = c(extreme = 1),
    scr_pthr = 0.05
  )

  expect_length(screened$betahat_orth, 1L)
  expect_true(is.finite(screened$betahat_orth))
  expect_true(is.finite(screened$sigma_orth))
})

test_that("correlation block threshold accepts the documented zero boundary", {
  traits <- c("trait_a", "trait_b")
  blocks <- fastASSETcli:::create_blocks(
    identity_correlation(traits), cor_thr = 0
  )
  expect_identical(blocks, list(traits))
})

test_that("published analysis standardization agrees on the multi-trait path", {
  skip_if_not_installed("ASSET")
  traits <- c("trait_positive", "trait_negative")
  beta <- c(0.35, -0.3)
  sigma <- c(0.1, 0.1)
  neff <- c(100, 400)
  correlation <- identity_correlation(traits)

  ours <- fastASSETcli:::fast_asset_2(
    snp = "rs_published_path",
    traits.lab = traits,
    beta.hat = beta,
    sigma.hat = sigma,
    Neff = neff,
    cor = correlation,
    block = list(),
    include_meta = TRUE
  )

  beta_published <- (beta / sigma) / sqrt(neff)
  sigma_published <- 1 / sqrt(neff)
  names(beta_published) <- names(sigma_published) <- traits
  screened <- fastASSETcli:::scr_transform(
    betahat = beta_published,
    SE = sigma_published,
    SNP = "rs_published_path",
    traits_i = traits,
    cor = correlation,
    block = list(),
    Neff = stats::setNames(neff, traits),
    scr_pthr = 0.05
  )
  reference <- ASSET::h.traits(
    snp.vars = "rs_published_path",
    traits.lab = names(screened$betahat_orth),
    beta.hat = screened$betahat_orth,
    sigma.hat = screened$sigma_orth,
    ncase = 2 / screened$sigma_orth^2,
    ncntl = 2 / screened$sigma_orth^2,
    cor = correlation[
      names(screened$betahat_orth), names(screened$betahat_orth), drop = FALSE
    ],
    side = 2,
    search = 2,
    cor.numr = FALSE,
    meta = TRUE,
    meth.pval = "DLM"
  )

  expect_identical(ours$qc$status, "PASS")
  expect_equal(
    ours$result$Subset.2sided$pval,
    reference$Subset.2sided$pval,
    tolerance = 1e-12
  )
  expect_equal(
    ours$result$Subset.2sided$pheno.1,
    reference$Subset.2sided$pheno.1
  )
  expect_equal(
    ours$result$Subset.2sided$pheno.2,
    reference$Subset.2sided$pheno.2
  )
})

test_that("Meta calculation does not alter directional subset results", {
  skip_if_not_installed("ASSET")
  traits <- c("trait_positive", "trait_negative")
  arguments <- list(
    snp = "rs_meta_invariance",
    traits.lab = traits,
    beta.hat = c(0.35, -0.3),
    sigma.hat = c(0.1, 0.1),
    Neff = c(100, 400),
    cor = identity_correlation(traits),
    block = list()
  )
  with_meta <- do.call(
    fastASSETcli:::fast_asset_2, c(arguments, list(include_meta = TRUE))
  )
  without_meta <- do.call(
    fastASSETcli:::fast_asset_2, c(arguments, list(include_meta = FALSE))
  )

  expect_equal(
    with_meta$result$Subset.2sided,
    without_meta$result$Subset.2sided,
    tolerance = 1e-12
  )
  expect_false(is.null(with_meta$result$Meta))
  expect_null(without_meta$result$Meta)
})

test_that("the ASSET fast_asset copy remains a missing-standardization control", {
  skip_if_not_installed("ASSET")
  asset_fast <- getExportedValue("ASSET", "fast_asset")
  body_text <- paste(deparse(body(asset_fast)), collapse = "\n")
  skip_if(
    grepl("sqrt\\(Neff\\)", body_text),
    "Installed ASSET copy now contains Neff standardization"
  )

  traits <- c("trait_a", "trait_b")
  correlation <- identity_correlation(traits)
  beta <- c(0.3, 0.4)
  sigma <- c(0.1, 0.1)
  neff <- c(100, 10000)
  ours <- fastASSETcli:::fast_asset_2(
    snp = "rs_negative_control",
    traits.lab = traits,
    beta.hat = beta,
    sigma.hat = sigma,
    Neff = neff,
    cor = correlation,
    block = list()
  )
  outlier <- asset_fast(
    snp = "rs_negative_control",
    traits.lab = traits,
    beta.hat = beta,
    sigma.hat = sigma,
    Neff = neff,
    cor = correlation,
    block = list(),
    scr_pthr = 0.05
  )

  expect_false(isTRUE(all.equal(
    as.numeric(ours$result$sigma.hat),
    as.numeric(outlier$sigma.hat),
    tolerance = 1e-12
  )))
})

test_that("pinned gqi fastASSET agrees on its packaged multi-trait fixture", {
  skip_if_not_installed("ASSET")
  skip_if_not_installed("fastASSET")
  pinned_sha <- "dbea4a67dd8699e1e1e7ca0e8d6277d9c602ffd2"
  description <- utils::packageDescription("fastASSET")
  installed_sha <- description[["RemoteSha"]]
  skip_if(
    is.null(installed_sha) || !identical(installed_sha, pinned_sha),
    paste0("fastASSET is not installed from pinned commit ", pinned_sha)
  )

  fixture <- new.env(parent = emptyenv())
  utils::data(
    "example_rs6678982", package = "fastASSET", envir = fixture
  )
  correlation <- fastASSETcli:::normalize_ldsc_intercept(
    fixture$ldscintmat, traits = fixture$traits
  )
  blocks <- fastASSETcli:::create_blocks(correlation)
  ours <- fastASSETcli:::fast_asset_2(
    snp = fixture$SNP,
    traits.lab = fixture$traits,
    beta.hat = fixture$betahat,
    sigma.hat = fixture$SE,
    Neff = fixture$Neff,
    cor = correlation,
    block = blocks,
    scr_pthr = 0.05,
    max_numtraits_per_side = 16L,
    include_meta = FALSE
  )
  upstream_fast <- getExportedValue("fastASSET", "fast_asset")
  upstream <- upstream_fast(
    snp = fixture$SNP,
    traits.lab = fixture$traits,
    beta.hat = fixture$betahat,
    sigma.hat = fixture$SE,
    Neff = fixture$Neff,
    cor = correlation,
    block = blocks,
    scr_pthr = 0.05,
    max_numtraits_per_side = 16L
  )

  expect_equal(
    ours$result$Subset.2sided$pval,
    upstream$Subset.2sided$pval,
    tolerance = 1e-12
  )
  expect_equal(
    ours$result$Subset.2sided$pheno.1,
    upstream$Subset.2sided$pheno.1
  )
  expect_equal(
    ours$result$Subset.2sided$pheno.2,
    upstream$Subset.2sided$pheno.2
  )
})
