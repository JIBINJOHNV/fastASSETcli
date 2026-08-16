# Reference-aligned fastASSET core.
#
# Important choices:
#   * The internal NEF field is the upstream fastASSET sample-size quantity:
#     total analyzed N for quantitative traits and Ncase*Ncontrol/(Ncase+
#     Ncontrol) for binary traits.
#   * ASSET is called with side = 2 and search = 2. Fixed-effect Meta output is
#     included by default but does not change the directional subset search.
#   * The default computational guard is 16 screened traits per direction.

assert_scalar_number <- function(x, name, lower = -Inf, upper = Inf,
                                 lower_open = FALSE, upper_open = FALSE) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(name, " must be one finite number.")
  }

  lower_bad <- if (lower_open) x <= lower else x < lower
  upper_bad <- if (upper_open) x >= upper else x > upper

  if (lower_bad || upper_bad) {
    left <- if (lower_open) "(" else "["
    right <- if (upper_open) ")" else "]"
    stop(name, " must be in ", left, lower, ", ", upper, right, ".")
  }

  invisible(TRUE)
}

normalize_ldsc_intercept <- function(intercept_matrix, traits = NULL,
                                     eigen_tolerance = 1e-8) {
  intercept_matrix <- as.matrix(intercept_matrix)
  storage.mode(intercept_matrix) <- "double"

  if (nrow(intercept_matrix) != ncol(intercept_matrix)) {
    stop("The LDSC intercept matrix must be square.")
  }
  if (nrow(intercept_matrix) == 0L) {
    stop("The LDSC intercept matrix is empty.")
  }
  if (is.null(rownames(intercept_matrix)) ||
      is.null(colnames(intercept_matrix))) {
    stop("The LDSC intercept matrix must have row and column trait names.")
  }
  if (anyNA(rownames(intercept_matrix)) || anyNA(colnames(intercept_matrix)) ||
      any(!nzchar(rownames(intercept_matrix))) ||
      any(!nzchar(colnames(intercept_matrix))) ||
      anyDuplicated(rownames(intercept_matrix)) ||
      anyDuplicated(colnames(intercept_matrix))) {
    stop("LDSC intercept row and column trait names must be non-empty and unique.")
  }
  if (!setequal(rownames(intercept_matrix), colnames(intercept_matrix))) {
    stop("The LDSC intercept matrix row and column names do not match.")
  }

  if (!is.null(traits)) {
    traits <- as.character(traits)
    if (anyNA(traits) || any(!nzchar(traits)) || anyDuplicated(traits)) {
      stop("Requested trait names must be non-empty and unique.")
    }
    missing_traits <- setdiff(traits, rownames(intercept_matrix))
    if (length(missing_traits) > 0L) {
      stop(
        "Traits missing from the LDSC intercept matrix: ",
        paste(missing_traits, collapse = ", ")
      )
    }
    intercept_matrix <- intercept_matrix[traits, traits, drop = FALSE]
  }

  if (any(!is.finite(intercept_matrix))) {
    stop("The LDSC intercept matrix contains non-finite values.")
  }

  intercept_matrix <- (intercept_matrix + t(intercept_matrix)) / 2
  intercept_variance <- diag(intercept_matrix)

  if (any(!is.finite(intercept_variance)) || any(intercept_variance <= 0)) {
    stop("Every diagonal LDSC intercept must be finite and greater than zero.")
  }

  # Trait-covariance normalization: D^(-1/2) I D^(-1/2).
  scale_factor <- sqrt(outer(intercept_variance, intercept_variance))
  correlation <- intercept_matrix / scale_factor
  correlation <- (correlation + t(correlation)) / 2
  diag(correlation) <- 1

  if (any(abs(correlation) > 1 + 1e-8)) {
    stop(
      "Normalizing the LDSC intercept matrix produced |correlation| > 1. ",
      "Inspect the LDSC intercept estimates; they are not a valid covariance matrix."
    )
  }

  # Remove only floating-point overshoot, never substantive invalid estimates.
  correlation[correlation > 1] <- 1
  correlation[correlation < -1] <- -1

  eigenvalues <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
  minimum_eigenvalue <- min(eigenvalues)
  if (!is.finite(minimum_eigenvalue) || minimum_eigenvalue <= eigen_tolerance) {
    stop(
      "The normalized trait correlation matrix is not positive definite. ",
      "Minimum eigenvalue = ", format(minimum_eigenvalue, digits = 8),
      ". Do not silently replace it with nearPD; inspect the LDSC inputs."
    )
  }

  attr(correlation, "minimum_eigenvalue") <- minimum_eigenvalue
  attr(correlation, "intercept_diagonal") <- intercept_variance
  correlation
}

scr_transform_in_block <- function(snp, traits.lab, beta.hat, sigma.hat,
                                   cor = NULL, n.eff = NULL, p.bound = 1) {
  k <- length(traits.lab)
  if (!is.numeric(beta.hat) || length(beta.hat) != k) {
    stop("beta.hat must be a numeric vector with one value per trait.")
  }
  if (!is.numeric(sigma.hat) || length(sigma.hat) != k) {
    stop("sigma.hat must be a numeric vector with one value per trait.")
  }

  names(beta.hat) <- traits.lab
  names(sigma.hat) <- traits.lab

  ord <- if (is.null(n.eff)) seq_len(k) else order(n.eff)
  beta.mat <- beta.hat[ord]
  sigma.mat <- sigma.hat[ord]
  t.lab <- traits.lab[ord]

  if (p.bound >= 1) {
    return(list(beta.scr = beta.mat, sigma.scr = sigma.mat))
  }

  cutoff <- stats::qnorm(p.bound / 2, lower.tail = FALSE)

  if (is.null(cor)) {
    U <- diag(k)
    Ui <- diag(k)
  } else {
    cor.ord <- cor[ord, ord, drop = FALSE]
    U <- chol(cor.ord)
    Ui <- solve(U)
  }

  z <- beta.mat / sigma.mat
  z.orthogonal <- as.vector(z %*% Ui)
  selected <- abs(z.orthogonal) > cutoff

  conditional_p <- 2 * stats::pnorm(
    abs(z.orthogonal), lower.tail = FALSE
  ) / p.bound
  conditional_p <- pmax(
    pmin(conditional_p, 0.9999999),
    .Machine$double.xmin
  )

  z.adjusted <- sign(z.orthogonal) * stats::qnorm(
    conditional_p / 2,
    lower.tail = FALSE
  )
  z.recorrelated <- as.vector(z.adjusted %*% U)

  beta.adjusted <- z.recorrelated * sigma.mat
  names(beta.adjusted) <- t.lab
  names(sigma.mat) <- t.lab
  names(selected) <- t.lab

  list(
    beta.scr = beta.adjusted[selected],
    sigma.scr = sigma.mat[selected]
  )
}

scr_transform <- function(betahat, SE, SNP, traits_i, cor, block, Neff,
                          scr_pthr = 0.05) {
  assert_scalar_number(
    scr_pthr, "scr_pthr", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )

  traits_i_block <- lapply(block, function(x) traits_i[traits_i %in% x])
  traits_i_block <- traits_i_block[lengths(traits_i_block) >= 2L]

  beta.adjusted <- numeric()
  sigma.adjusted <- numeric()

  if (length(traits_i_block) > 0L) {
    for (block_traits in traits_i_block) {
      screened <- scr_transform_in_block(
        snp = SNP,
        traits.lab = block_traits,
        beta.hat = betahat[block_traits],
        sigma.hat = SE[block_traits],
        n.eff = Neff[block_traits],
        cor = cor[block_traits, block_traits, drop = FALSE],
        p.bound = scr_pthr
      )
      beta.adjusted <- c(beta.adjusted, screened$beta.scr)
      sigma.adjusted <- c(sigma.adjusted, screened$sigma.scr)
    }
  }

  blocked_traits <- unique(unlist(traits_i_block, use.names = FALSE))
  independent_traits <- setdiff(traits_i, blocked_traits)

  if (length(independent_traits) > 0L) {
    beta.independent <- betahat[independent_traits]
    se.independent <- SE[independent_traits]
    raw_p <- 2 * stats::pnorm(
      abs(beta.independent / se.independent),
      lower.tail = FALSE
    )
    selected <- raw_p < scr_pthr

    if (any(selected)) {
      beta.selected <- beta.independent[selected]
      se.selected <- se.independent[selected]
      conditional_p <- raw_p[selected] / scr_pthr
      conditional_p <- pmax(
        pmin(conditional_p, 0.9999999),
        .Machine$double.xmin
      )
      beta.selected <- se.selected * sign(beta.selected) * stats::qnorm(
        conditional_p / 2,
        lower.tail = FALSE
      )

      beta.adjusted <- c(beta.adjusted, beta.selected)
      sigma.adjusted <- c(sigma.adjusted, se.selected)
    }
  }

  list(
    betahat_orth = beta.adjusted,
    sigma_orth = sigma.adjusted
  )
}

make_fastasset_qc <- function(snp, status, severity, message,
                              n_input_traits, valid_traits, input_traits = NULL,
                              screened_positive = character(),
                              screened_negative = character(),
                              screening_applied = NA,
                              analysis_scope = "NOT_ANALYZED") {
  n_positive <- length(screened_positive)
  n_negative <- length(screened_negative)
  invalid_traits <- if (is.null(input_traits)) {
    character()
  } else {
    setdiff(input_traits, valid_traits)
  }

  list(
    ID = as.character(snp),
    status = status,
    severity = severity,
    message = message,
    screening_applied = as.logical(screening_applied),
    analysis_scope = as.character(analysis_scope),
    n_input_traits = as.integer(n_input_traits),
    n_valid_traits = as.integer(length(valid_traits)),
    n_invalid_traits = as.integer(n_input_traits - length(valid_traits)),
    invalid_traits = paste(invalid_traits, collapse = ";"),
    n_screened_traits = as.integer(n_positive + n_negative),
    n_screened_positive = as.integer(n_positive),
    n_screened_negative = as.integer(n_negative),
    screened_positive_traits = paste(screened_positive, collapse = ";"),
    screened_negative_traits = paste(screened_negative, collapse = ";"),
    candidate_subsets_positive = 2^n_positive - 1,
    candidate_subsets_negative = 2^n_negative - 1,
    candidate_subsets_total = (2^n_positive - 1) + (2^n_negative - 1)
  )
}

make_single_trait_fastasset_result <- function(snp, trait, beta, sigma,
                                                include_meta = TRUE) {
  if (length(trait) != 1L || length(beta) != 1L || length(sigma) != 1L ||
      !is.finite(beta) || !is.finite(sigma) || sigma <= 0) {
    stop("Single-trait closed-form input must contain one finite effect and positive SE.")
  }

  snp <- as.character(snp)
  trait <- as.character(trait)
  beta <- as.numeric(beta)
  sigma <- as.numeric(sigma)
  p_two_sided <- as.numeric(2 * stats::pnorm(-abs(beta / sigma)))
  positive <- beta >= 0

  pheno_positive <- matrix(
    positive,
    nrow = 1L,
    dimnames = list(snp, trait)
  )
  pheno_negative <- matrix(
    !positive,
    nrow = 1L,
    dimnames = list(snp, trait)
  )

  positive_value <- function(value) if (positive) value else NA_real_
  negative_value <- function(value) if (positive) NA_real_ else value

  subset_two_sided <- list(
    pval = stats::setNames(p_two_sided, snp),
    pval.1 = stats::setNames(if (positive) p_two_sided else 1, snp),
    beta.1 = stats::setNames(positive_value(beta), snp),
    sd.1 = stats::setNames(positive_value(sigma), snp),
    pheno.1 = pheno_positive,
    pval.2 = stats::setNames(if (positive) 1 else p_two_sided, snp),
    beta.2 = stats::setNames(negative_value(beta), snp),
    sd.2 = stats::setNames(negative_value(sigma), snp),
    pheno.2 = pheno_negative,
    sd.1.meta = stats::setNames(positive_value(sigma), snp),
    sd.2.meta = stats::setNames(negative_value(sigma), snp)
  )

  meta <- if (include_meta) {
    list(
      pval = stats::setNames(p_two_sided, snp),
      beta = stats::setNames(beta, snp),
      sd = stats::setNames(sigma, snp)
    )
  } else {
    NULL
  }

  list(
    Meta = meta,
    Subset.1sided = NULL,
    Subset.2sided = subset_two_sided,
    snp.vars = snp,
    traits.lab = trait,
    beta.hat = stats::setNames(beta, trait),
    sigma.hat = stats::setNames(sigma, trait),
    search = 2,
    side = 2,
    meta = include_meta,
    which = "single_trait_closed_form"
  )
}

direction_limit_exceeded <- function(beta, max_numtraits_per_side) {
  sum(beta >= 0) > max_numtraits_per_side ||
    sum(beta < 0) > max_numtraits_per_side
}

# Run fastASSET for one SNP and return both the ASSET object and auditable QC.
fast_asset_2 <- function(snp, traits.lab, beta.hat, sigma.hat, Neff, cor, block,
                         scr_pthr = 0.05, max_numtraits_per_side = 16L,
                         min_available_traits = 1L, meth_pval = "DLM",
                         include_meta = TRUE) {
  n_traits <- length(traits.lab)
  if (any(c(length(beta.hat), length(sigma.hat), length(Neff)) != n_traits)) {
    stop("traits.lab, beta.hat, sigma.hat and Neff must have equal lengths.")
  }
  if (anyDuplicated(traits.lab)) {
    stop("traits.lab contains duplicate trait names.")
  }

  assert_scalar_number(
    scr_pthr, "scr_pthr", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  assert_scalar_number(max_numtraits_per_side, "max_numtraits_per_side", 1, 30)
  assert_scalar_number(min_available_traits, "min_available_traits", 1, n_traits)
  if (length(include_meta) != 1L || is.na(include_meta) || !is.logical(include_meta)) {
    stop("include_meta must be TRUE or FALSE.")
  }

  cor <- as.matrix(cor)
  if (is.null(rownames(cor)) || is.null(colnames(cor))) {
    if (!all(dim(cor) == c(n_traits, n_traits))) {
      stop("Unnamed cor must have one row and column per input trait.")
    }
    rownames(cor) <- colnames(cor) <- traits.lab
  }
  missing_cor_traits <- setdiff(traits.lab, rownames(cor))
  if (length(missing_cor_traits) > 0L) {
    stop(
      "Traits missing from cor: ",
      paste(missing_cor_traits, collapse = ", ")
    )
  }

  valid <- is.finite(beta.hat) & is.finite(sigma.hat) & is.finite(Neff) &
    sigma.hat > 0 & Neff > 0
  valid_traits <- traits.lab[valid]

  if (length(valid_traits) == 0L) {
    qc <- make_fastasset_qc(
      snp, "NO_VALID_TRAITS", "WARNING",
      "No trait had finite BETA, SE and NEF with SE > 0 and NEF > 0.",
      n_traits, valid_traits, traits.lab,
      screening_applied = FALSE,
      analysis_scope = "NOT_ANALYZED"
    )
    return(list(result = NULL, qc = qc))
  }

  if (length(valid_traits) < min_available_traits) {
    qc <- make_fastasset_qc(
      snp, "INSUFFICIENT_VALID_TRAITS", "WARNING",
      paste0(
        "Only ", length(valid_traits), " valid trait(s); minimum is ",
        min_available_traits, "."
      ),
      n_traits, valid_traits, traits.lab,
      screening_applied = FALSE,
      analysis_scope = "NOT_ANALYZED"
    )
    return(list(result = NULL, qc = qc))
  }

  beta.hat <- beta.hat[valid]
  sigma.hat <- sigma.hat[valid]
  Neff <- Neff[valid]
  names(beta.hat) <- names(sigma.hat) <- names(Neff) <- valid_traits

  # Reference transformation. Neff follows the upstream fastASSET definition.
  beta.standardized <- beta.hat / sigma.hat / sqrt(Neff)
  sigma.standardized <- 1 / sqrt(Neff)

  # The published genome-wide analysis bypasses screening when only one valid
  # trait exists and uses the ordinary two-sided GWAS z test.  The standardized
  # beta/SE ratio is algebraically identical to the original beta/SE ratio.
  if (length(valid_traits) == 1L) {
    result <- make_single_trait_fastasset_result(
      snp = snp,
      trait = valid_traits,
      beta = beta.standardized,
      sigma = sigma.standardized,
      include_meta = include_meta
    )
    qc <- make_fastasset_qc(
      snp, "SINGLE_VALID_TRAIT_NO_SCREEN", "OK",
      paste0(
        "Only trait ", valid_traits,
        " was valid; pre-screening and ASSET enumeration were bypassed and ",
        "the published closed-form two-sided z test was used."
      ),
      n_traits, valid_traits, traits.lab,
      screening_applied = FALSE,
      analysis_scope = "SINGLE_VALID_TRAIT_UNSCREENED"
    )
    return(list(result = result, qc = qc))
  }

  cormat <- cor[valid_traits, valid_traits, drop = FALSE]
  screened <- scr_transform(
    betahat = beta.standardized,
    SE = sigma.standardized,
    SNP = snp,
    traits_i = valid_traits,
    cor = cormat,
    block = block,
    Neff = Neff,
    scr_pthr = scr_pthr
  )

  screened_beta <- screened$betahat_orth
  screened_se <- screened$sigma_orth
  screened_traits <- names(screened_beta)

  # This exactly follows the split actually used by ASSET::h.traits2:
  # non-negative adjusted effects versus negative adjusted effects.
  positive_traits <- screened_traits[screened_beta >= 0]
  negative_traits <- screened_traits[screened_beta < 0]

  if (length(screened_traits) == 0L) {
    qc <- make_fastasset_qc(
      snp, "NO_TRAIT_PASSED_SCREEN", "INFO",
      paste0("No valid trait had pre-screening P < ", scr_pthr, "."),
      n_traits, valid_traits, traits.lab,
      screening_applied = TRUE,
      analysis_scope = "NOT_ANALYZED"
    )
    return(list(result = NULL, qc = qc))
  }

  if (direction_limit_exceeded(screened_beta, max_numtraits_per_side)) {
    qc <- make_fastasset_qc(
      snp, "DIRECTION_LIMIT_EXCEEDED", "HIGH",
      paste0(
        "Screened counts were positive=", length(positive_traits),
        " and negative=", length(negative_traits),
        "; configured maximum is ", max_numtraits_per_side,
        " per direction. SNP was not sent to exhaustive subset search."
      ),
      n_traits, valid_traits, traits.lab, positive_traits, negative_traits,
      screening_applied = TRUE,
      analysis_scope = "SCREENED_DIRECTION_LIMIT_EXCLUDED"
    )
    return(list(result = NULL, qc = qc))
  }

  # The published analysis uses the adjusted one-trait z statistic directly
  # instead of routing it through ASSET's DLM machinery.
  if (length(screened_traits) == 1L) {
    result <- make_single_trait_fastasset_result(
      snp = snp,
      trait = screened_traits,
      beta = screened_beta,
      sigma = screened_se,
      include_meta = include_meta
    )
    qc <- make_fastasset_qc(
      snp, "SINGLE_TRAIT_AFTER_SCREEN", "OK",
      paste0(
        "Only trait ", screened_traits,
        " remained after pre-screening; the published closed-form two-sided ",
        "test was applied to its selection-adjusted z statistic."
      ),
      n_traits, valid_traits, traits.lab, positive_traits, negative_traits,
      screening_applied = TRUE,
      analysis_scope = "SINGLE_SCREENED_SELECTION_ADJUSTED_TRAIT"
    )
    return(list(result = result, qc = qc))
  }

  htraits_args <- list(
    snp.vars = as.character(snp),
    traits.lab = screened_traits,
    beta.hat = screened_beta,
    sigma.hat = screened_se,
    ncase = 2 / screened_se^2,
    ncntl = 2 / screened_se^2,
    side = 2,
    search = 2,
    meta = include_meta,
    meth.pval = meth_pval
  )

  htraits_args$cor <- cormat[screened_traits, screened_traits, drop = FALSE]
  htraits_args$cor.numr <- FALSE

  if (!requireNamespace("ASSET", quietly = TRUE)) {
    stop("The ASSET package is required for multi-trait subset enumeration.")
  }
  result <- do.call(ASSET::h.traits, htraits_args)
  qc <- make_fastasset_qc(
    snp, "PASS", "OK", "fastASSET completed.",
    n_traits, valid_traits, traits.lab, positive_traits, negative_traits,
    screening_applied = TRUE,
    analysis_scope = "SCREENED_SELECTION_ADJUSTED_TRAITS"
  )

  list(result = result, qc = qc)
}

extract_fastasset_meta <- function(fastasset_outcome) {
  result <- fastasset_outcome$result
  qc <- fastasset_outcome$qc

  if (is.null(result) || is.null(result$Meta)) return(NULL)

  meta_result <- result$Meta
  analyzed_traits <- as.character(result$traits.lab)
  scalar <- function(x) {
    if (is.null(x) || length(x) == 0L) NA_real_ else as.numeric(x[1L])
  }

  list(
    ID = qc$ID,
    status = qc$status,
    severity = qc$severity,
    screening_applied = qc$screening_applied,
    meta_scope = qc$analysis_scope,
    n_analyzed_traits = as.integer(length(analyzed_traits)),
    analyzed_traits = paste(analyzed_traits, collapse = ";"),
    p_meta = scalar(meta_result$pval),
    beta_meta_adjusted = scalar(meta_result$beta),
    se_meta_adjusted = scalar(meta_result$sd),
    n_screened_traits = qc$n_screened_traits,
    n_screened_positive = qc$n_screened_positive,
    n_screened_negative = qc$n_screened_negative,
    screened_positive_traits = qc$screened_positive_traits,
    screened_negative_traits = qc$screened_negative_traits
  )
}

extract_fastasset_two_sided <- function(fastasset_outcome) {
  result <- fastasset_outcome$result
  qc <- fastasset_outcome$qc

  if (is.null(result) || is.null(result$Subset.2sided)) {
    return(NULL)
  }

  subset_result <- result$Subset.2sided

  scalar <- function(x) {
    if (is.null(x) || length(x) == 0L) NA_real_ else as.numeric(x[1L])
  }

  selected_traits <- function(indicator) {
    if (is.null(indicator)) return(character())
    if (is.matrix(indicator)) {
      values <- as.logical(indicator[1L, , drop = TRUE])
      trait_names <- colnames(indicator)
    } else {
      values <- as.logical(indicator)
      trait_names <- names(indicator)
    }
    if (is.null(trait_names)) trait_names <- result$traits.lab
    trait_names[which(values)]
  }

  selected_positive <- selected_traits(subset_result$pheno.1)
  selected_negative <- selected_traits(subset_result$pheno.2)

  list(
    ID = qc$ID,
    status = qc$status,
    severity = qc$severity,
    screening_applied = qc$screening_applied,
    analysis_scope = qc$analysis_scope,
    p_two_sided = scalar(subset_result$pval),
    p_positive = scalar(subset_result$pval.1),
    p_negative = scalar(subset_result$pval.2),
    beta_positive_adjusted = scalar(subset_result$beta.1),
    se_positive_adjusted = scalar(subset_result$sd.1),
    se_positive_meta = scalar(subset_result$sd.1.meta),
    beta_negative_adjusted = scalar(subset_result$beta.2),
    se_negative_adjusted = scalar(subset_result$sd.2),
    se_negative_meta = scalar(subset_result$sd.2.meta),
    n_selected_positive = length(selected_positive),
    n_selected_negative = length(selected_negative),
    selected_positive_traits = paste(selected_positive, collapse = ";"),
    selected_negative_traits = paste(selected_negative, collapse = ";"),
    n_screened_positive = qc$n_screened_positive,
    n_screened_negative = qc$n_screened_negative,
    screened_positive_traits = qc$screened_positive_traits,
    screened_negative_traits = qc$screened_negative_traits,
    candidate_subsets_total = qc$candidate_subsets_total
  )
}

create_blocks <- function(cormat, cor_thr = 0.2) {
  cormat <- as.matrix(cormat)
  assert_scalar_number(
    cor_thr, "cor_thr", lower = 0, upper = 1,
    lower_open = FALSE, upper_open = TRUE
  )

  if (nrow(cormat) != ncol(cormat)) stop("cormat must be square.")
  if (nrow(cormat) < 2L) return(list())
  if (is.null(rownames(cormat))) stop("cormat must have trait names.")

  correlation_distance <- stats::as.dist(1 - abs(cormat))
  hc <- stats::hclust(correlation_distance)
  tree <- stats::cutree(hc, h = 1 - cor_thr)
  retained_groups <- as.integer(names(table(tree))[table(tree) >= 2L])
  lapply(retained_groups, function(group) names(tree)[tree == group])
}
