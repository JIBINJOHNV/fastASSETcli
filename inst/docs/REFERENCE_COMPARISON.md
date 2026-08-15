# Reference comparison and corrections

This review compares the supplied `fastASSET.R`, `fast_asset_analysis.R`, and
`helpers.R` with the upstream fastASSET implementation, the ASSET source/manual,
the GenomicSEM LDSC object returned by `ldsc()`, and the fastASSET publication.

## Corrections made

| Area | Supplied code | Corrected code | Severity | If left unchanged |
|---|---|---|---|---|
| Direction search | Called `h.traits(..., search=NULL)`. | Calls `side=2, search=2`, matching upstream fastASSET. | **CRITICAL** | ASSET also performs the one-sided all-subset search. Runtime rises sharply and outputs no longer represent only the intended positive/negative split search. |
| One-direction guard | Default was changed from 16 to 100. | Default restored to 16 screened traits per direction; CLI permits 1–30 and labels values above 16 as sensitivity settings. | **CRITICAL** | A SNP with 100 screened effects in one direction implies `2^100-1` candidate subsets and is not computationally feasible, regardless of 90 CPUs or 900 GB RAM. |
| Screening threshold | The function accepted `scr_pthr` but passed a hard-coded `0.05` into `scr_transform()`. | The CLI value is passed through exactly and written to settings. | **HIGH** | A user-requested threshold would silently be ignored, invalidating sensitivity analyses and reproducibility. |
| Direction definition | The supplied guard used positive `>0` and negative `<=0`. | Uses non-negative `>=0` and negative `<0`, exactly as `ASSET::h.traits2`. | **MEDIUM** | A zero adjusted effect could be counted on the opposite side from the actual ASSET search. Usually rare, but QC and guard decisions could disagree with ASSET. |
| Missing/invalid values | Checked only non-missing Beta and SE before transformation. | Per SNP, requires finite Beta/SE/NEF, `SE>0`, and `NEF>0`; invalid traits are recorded without dropping the entire SNP. | **HIGH** | Infinite values, zero SE, or invalid NEF can produce invalid transformed statistics, crashes, or silent non-finite results. |
| Quantitative sample size | The supplied fastASSET function already used `Neff` directly, but the wide-table contract did not clearly enforce quantitative NEF columns and valid values. | `<trait>.NEF` is explicitly required and used directly as total analyzed quantitative sample size; no case/control conversion is applied. | **INFO / REQUIREMENT** | Ambiguous sample-size columns make it easier to supply the wrong quantity; an incorrect value changes the transformed beta/SE, screening, subset membership, and P values. |
| LDSC covariance scale | Forced `diag(I)=1` and used off-diagonal intercepts unchanged. | Computes `D^(-1/2) I D^(-1/2)`, symmetrizes, checks bounds, and requires positive definiteness. | **CRITICAL** | Off-diagonal covariance estimates remain on the wrong scale when univariate intercepts are not one, so the correlation supplied to screening and ASSET is incorrect. |
| LDSC trait order | The supplied pipeline generated LDSC internally and copied names from `S`; it did not provide the requested audited loader for a separate saved GenomicSEM object. | Existing-object mode extracts `$I`, recovers labels from matching `I`, `S`, or `S_Stand` dimnames, and reorders both axes. Automatic mode first requires the manifest trait set to equal the `.Beta` trait set, reorders the manifest before munging/LDSC, assigns both-axis dimnames to generated `I` and `S`, saves the object separately, and then repeats the extraction audit. | **CRITICAL** | If a supplied or generated matrix is used without verified two-axis reordering, correlations can be attached to the wrong traits, producing scientifically invalid results without necessarily causing an R error. |
| Automatic LDSC resume | The supplied code reused an LDSC `.RData` whenever the expected path existed. | Generated output is stored in a signature-addressed directory based on trait order, files, reference metadata, prevalence values and LDSC settings; reuse additionally requires a completion marker. | **HIGH** | A partially written or scientifically incompatible LDSC object could be silently reused. |
| Positive definiteness | Did not fully protect the matrix transformation path. | Rejects non-finite, out-of-range, or non-positive-definite normalized matrices; never silently calls `nearPD`. | **HIGH** | Cholesky decomposition can fail inside screening, or a silently altered matrix can change the analysis without a transparent scientific decision. |
| Meta-analysis | Meta output was formatted through `h.summary()` and Beta/SE could be reconstructed from rounded OR confidence intervals. | Calls `h.traits(..., meta=TRUE)` and extracts raw `Meta$pval`, `Meta$beta`, and `Meta$sd` directly. | **HIGH** | Reconstructing estimates from rounded display values loses precision and may not reproduce the actual ASSET result. |
| Meta scope | The meaning of the meta result was not explicit. | The meta file states `SCREENED_SELECTION_ADJUSTED_TRAITS` and records all screened trait names/counts. | **HIGH** | Users may incorrectly interpret the result as an ordinary fixed-effect meta-analysis of all 200–250 original GWAS effects. |
| Full standard ASSET stage | Attempted ordinary `h.traits` across the complete trait set. | Removed the exhaustive all-200–250-trait subset stage. | **CRITICAL** | Exhaustive subset enumeration is computationally intractable; more CPUs parallelize SNPs, not the subset search within one SNP. |
| Resume integrity | Output existence checks could reuse partially incompatible files. | A run signature includes input/LDSC file metadata, trait order, correlation matrix, and scientific settings; a chunk resumes only when result, meta, QC, and completion files all exist. | **HIGH** | Results from an earlier parameter or matrix version can be mixed into a later run. |
| Error visibility | Several SNP/chunk failures could be difficult to distinguish from ordinary null results. | Every input SNP receives status, severity, message, directional counts, candidate counts, and runtime. | **HIGH** | Missing results may be mistaken for non-significant results, biasing downstream completeness and interpretation. |

## What “Meta” means in this bundle

For each SNP, fastASSET first standardizes the quantitative inputs using NEF,
pre-screens traits, and adjusts the retained statistics for that selection.
ASSET then receives only those retained adjusted traits. With `meta=TRUE`, its
`Meta` component is the fixed-effect meta-analysis of that screened,
selection-adjusted set. It is not a meta-analysis of every original trait.

This is kept separate from the two-direction subset result:

- `*_fastasset_results.tsv.gz` contains the positive/negative subset search.
- `*_fastasset_meta.tsv.gz` contains ASSET's fixed-effect `Meta` component.

## Primary settings for this analysis

- Traits at input: 200–250.
- Quantitative NEF: total analyzed sample size, used directly.
- CPUs: `--ncores 90` on Linux/macOS.
- Primary screening threshold: `--scr-pthr 0.05`.
- Primary direction guard: `--max-traits-per-side 16`.
- Meta output: `--include-meta TRUE` (also the default).

## Reference behavior used

- Upstream fastASSET standardizes as `beta / SE / sqrt(Neff)`, uses
  `SE = 1 / sqrt(Neff)`, pre-screens, and calls `h.traits` with `side=2` and
  `search=2`.
- ASSET `h.traits2` searches the non-negative and negative adjusted-effect
  groups separately.
- ASSET `meta=TRUE` returns `Meta = list(pval, beta, sd)`.
- GenomicSEM `ldsc()` returns a list containing `I` and `S`; trait names are
  assigned to `S`, so `S` is a necessary fallback when `I` has no usable names.

## Sources

- fastASSET source and tutorial: <https://github.com/gqi/fastASSET>
- ASSET source: <https://github.com/sbstatgen/ASSET>
- ASSET Bioconductor manual: <https://bioconductor.org/packages/release/bioc/manuals/ASSET/man/ASSET.pdf>
- GenomicSEM LDSC source: <https://github.com/GenomicSEM/GenomicSEM/blob/master/R/ldsc.R>
- Qi et al. 2024: <https://www.nature.com/articles/s41467-024-51075-5>
