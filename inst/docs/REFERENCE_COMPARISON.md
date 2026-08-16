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
| Single valid trait before screening | The wrapper could skip the SNP through a two-trait minimum or send one trait through screening. | The default minimum is one. A lone valid trait bypasses screening and ASSET enumeration and uses the published closed-form two-sided GWAS z test. | **HIGH** | A valid univariate result is lost, or a trait that the published workflow always analyzes is incorrectly removed by screening. |
| Single trait after screening | The wrapper routed every retained trait count through the same ASSET path. | A lone screened trait uses the published closed-form test on its selection-adjusted z statistic, with a dedicated QC status. | **HIGH** | A simple univariate result is unnecessarily routed through DLM machinery and its provenance cannot be distinguished from a multi-trait subset search. |
| Genome-wide direction-limit handling | The standalone package stops, while another ASSET copy only prints a warning. | Records `DIRECTION_LIMIT_EXCEEDED`, excludes the SNP from enumeration and continues, exactly as the released paper analysis script. | **HIGH** | Stopping loses the remainder of a genome-wide run; merely printing can allow an infeasible subset search to proceed. |
| Missing/invalid values | Checked only non-missing Beta and SE before transformation. | Per SNP, requires finite Beta/SE/NEF, `SE>0`, and `NEF>0`; invalid traits are recorded without dropping the entire SNP. | **HIGH** | Infinite values, zero SE, or invalid NEF can produce invalid transformed statistics, crashes, or silent non-finite results. |
| Analysis input contract | Required a separately prepared multi-trait wide table in addition to the per-trait manifest. | The manifest is the only summary-statistic input. Its per-trait TSV/VCF files are validated and combined internally into the exact wide FastASSET representation. | **HIGH** | Requiring two independently prepared representations permits trait-order, SNP, effect, SE, sample-size, and allele-orientation disagreements between LDSC and FastASSET. |
| Quantitative sample size | The supplied fastASSET function already used `Neff` directly, but the wide-table contract did not clearly enforce quantitative NEF columns and valid values. | `<trait>.NEF` is explicitly required and used directly as total analyzed quantitative sample size; no case/control conversion is applied. | **INFO / REQUIREMENT** | Ambiguous sample-size columns make it easier to supply the wrong quantity; an incorrect value changes the transformed beta/SE, screening, subset membership, and P values. |
| Binary fastASSET sample size | No single manifest-driven construction rule was implemented. | Uses the upstream fastASSET definition `Ncase*Ncontrol/(Ncase+Ncontrol)`. If counts are unavailable but total N and sample prevalence are provided, uses the algebraically equivalent `N*sample.prev*(1-sample.prev)`. | **CRITICAL** | Using total case-control N, or a four-times-rescaled effective N, changes the beta/SE transformation and therefore screening and subset P values. |
| Cross-trait effect orientation | A separately prepared wide table had no allele columns, so its effect alignment could not be audited. | For tabular inputs, the SNP-ID union is filled by indexed matching and per-trait A1/A2 are classified as exact, swapped or incompatible; swaps reverse BETA and incompatible observations become `NA` and are logged. For VCF inputs, the union key is ordered `CHR_POS_REF_ALT`, so allele identity is already part of the key and the redundant comparison/flip is skipped. | **CRITICAL** | A sign error can move a trait from the positive to the negative search and fundamentally change the selected subset and two-sided result. |
| GWAS-VCF effect orientation | No VCF input conversion was implemented. | For manifest VCF/BCF inputs, `FORMAT/ES` is explicitly interpreted relative to ALT, so the generated LDSC table uses `A1=ALT` and `A2=REF`. | **CRITICAL** | Swapping A1/A2 without reversing the effect changes its sign and can invert genetic covariance estimates. |
| GWAS-VCF identifier routing | No VCF branch was implemented. | Retains original VCF IDs, including repeated rsIDs, for HapMap3 ID-and-allele matching in `GenomicSEM::munge()`. Separately creates `FASTASSET_ID=CHR_POS_REF_ALT` for the multi-trait union and writes a complete rsID mapping file. | **HIGH / REQUIREMENT** | Rejecting duplicate rsIDs before munging can discard valid multiallelic records; merging FastASSET on rsID alone can collapse distinct alleles into one row. |
| GWAS-VCF P value | No conversion from GWAS-VCF `LP` was implemented. | Generates raw `P=10^(-LP)` and caps only values below the smallest positive representable R number. | **HIGH** | Passing `LP` as though it were raw P produces invalid Z statistics; allowing numeric underflow to zero produces infinite Z statistics. |
| Binary LDSC sample size and prevalence | Binary VCF counts were not extracted. | When `POPULATION_PREV` is supplied, uses per-SNP total `N=NC+NCO` and derives one trait-level `SAMPLE_PREV=median(NC)/(median(NC)+median(NCO))` across valid SNP rows. Header totals are retained only as provenance. | **CRITICAL** | Liability-scale conversion and LDSC regression weights are wrong when prevalence or sample size is wrong. |
| Quantitative VCF sample size | No VCF rule was implemented. | With both prevalences absent, uses `FORMAT/NEF` directly as the quantitative total sample size, as required for the supplied files. | **HIGH** | Treating a quantitative total N as a case-control effective N, or applying another conversion, changes LDSC weights and downstream covariance estimates. |
| LDSC covariance scale | Forced `diag(I)=1` and used off-diagonal intercepts unchanged. | Computes `D^(-1/2) I D^(-1/2)`, symmetrizes, checks bounds, and requires positive definiteness. This is the transformation used by the released 116-trait paper analysis before clustering, screening and ASSET. | **CRITICAL** | Off-diagonal covariance estimates remain on the wrong scale when univariate intercepts are not one, so the correlation supplied to screening and ASSET is incorrect. |
| LDSC/manifest trait matching | Native `GenomicSEM::ldsc()` constructs `$I` and `$S` in the same order and normally exposes the trait names on `$S`; the earlier loader did not implement the requested partial-intersection workflow. | Existing-object mode uses `$S` names to label both `$I` axes, compares them with the manifest, warns for traits missing from either source, retains only common traits in manifest order, and subsets the manifest rows—including VCF `SAMPLE`—before input preparation. A complete match audit is saved. Automatic mode remains exact manifest order. | **REQUIREMENT** | Without intersection and audited manifest-order subsetting, an unmatched manifest trait stops the run or a non-audited manual workaround can apply the wrong file/sample to a matrix position. |
| Automatic LDSC resume | The supplied code reused an LDSC `.RData` whenever the expected path existed. | Generated output is stored in a signature-addressed directory based on trait order, files, reference metadata, prevalence values and LDSC settings; reuse additionally requires a completion marker. | **HIGH** | A partially written or scientifically incompatible LDSC object could be silently reused. |
| Positive definiteness | Did not fully protect the matrix transformation path. | Rejects non-finite, out-of-range, or non-positive-definite normalized matrices; never silently calls `nearPD`. | **HIGH** | Cholesky decomposition can fail inside screening, or a silently altered matrix can change the analysis without a transparent scientific decision. |
| Meta-analysis | Meta output was formatted through `h.summary()` and Beta/SE could be reconstructed from rounded OR confidence intervals. | Calls `h.traits(..., meta=TRUE)` and extracts raw `Meta$pval`, `Meta$beta`, and `Meta$sd` directly. | **HIGH** | Reconstructing estimates from rounded display values loses precision and may not reproduce the actual ASSET result. |
| Meta scope | The meaning of the meta result was not explicit. | The meta file records the analysis scope: `SCREENED_SELECTION_ADJUSTED_TRAITS` for the multi-trait path, `SINGLE_SCREENED_SELECTION_ADJUSTED_TRAIT` for one trait remaining after screening, or `SINGLE_VALID_TRAIT_UNSCREENED` when only one trait was valid before screening. It also records the analyzed trait names/counts. | **HIGH** | Users may incorrectly interpret the result as an ordinary fixed-effect meta-analysis of all 200–250 original GWAS effects. |
| Full standard ASSET stage | Attempted ordinary `h.traits` across the complete trait set. | Removed the exhaustive all-200–250-trait subset stage. | **CRITICAL** | Exhaustive subset enumeration is computationally intractable; more CPUs parallelize SNPs, not the subset search within one SNP. |
| Resume integrity | Output existence checks could reuse partially incompatible files. | A run signature includes input/LDSC file metadata, trait order, correlation matrix, and scientific settings; a chunk resumes only when result, meta, QC, direction-limit report, and completion files all exist. | **HIGH** | Results from an earlier parameter or matrix version can be mixed into a later run. |
| Error visibility | Several SNP/chunk failures could be difficult to distinguish from ordinary null results. | Every input SNP receives status, severity, message, directional counts, candidate counts, and runtime. Direction-limit exclusions are additionally combined into a dedicated file while analysis continues for all other SNPs. | **HIGH** | Missing results may be mistaken for non-significant results, biasing downstream completeness and interpretation. |

## What “Meta” means in this bundle

For each SNP, fastASSET first standardizes the quantitative inputs using NEF,
pre-screens traits, and adjusts the retained statistics for that selection.
ASSET then receives only those retained adjusted traits. With `meta=TRUE`, its
`Meta` component is the fixed-effect meta-analysis of that screened,
selection-adjusted set. It is not a meta-analysis of every original trait.

This is kept separate from the two-direction subset result:

- `*_fastasset_results.tsv.gz` contains the positive/negative subset search.
- `*_fastasset_meta.tsv.gz` contains ASSET's fixed-effect `Meta` component for
  multi-trait rows and the equivalent closed-form result for one-trait rows.

For a single valid trait, there is no multi-trait subset enumeration. The Meta
file contains the equivalent one-trait result, and `meta_scope` distinguishes
an unscreened single valid trait from a single selection-adjusted screened
trait. `n_analyzed_traits` and `analyzed_traits` identify exactly which traits
contributed to every Meta row.

## Five-source provenance

The upstream artifacts are not behaviorally identical. “Reference-aligned” in
this package means following the released genome-wide analysis where it gives a
clear production rule, following `ASSET::h.traits2` for the actual direction
split, and retaining explicit safety corrections.

| Behavior | README tutorial | `fastASSET::fast_asset` | Released paper analysis | ASSET repository `fast_asset` | `ASSET::h.traits2` | fastASSETcli 0.5.1 |
|---|---|---|---|---|---|---|
| LDSC matrix scale | Passes the supplied example `ldscintmat` directly | Requires `cor` to be a z-statistic correlation with diagonal one | Explicitly applies `D^(-1/2) I D^(-1/2)` before every consumer | Requires correlation-scale `cor`; no normalization | Consumes the supplied correlation | Normalizes GenomicSEM covariance before every consumer |
| Neff standardization | Uses the package | `beta/SE/sqrt(Neff)` and `SE=1/sqrt(Neff)` | `z/sqrt(N)` and `SE=1/sqrt(N)`, algebraically identical | Missing | Expects prepared effects, SEs and sample-size inputs | Matches standalone package and paper script |
| One valid trait before screening | Uses the package | Calls `h.traits` without screening; the standalone source contains an uppercase `SNP` variable bug in this branch | Closed-form two-sided test on the original z statistic | Calls `h.traits` without screening | General one-trait calculation | Published closed form; `SINGLE_VALID_TRAIT_NO_SCREEN` |
| One trait after screening | Uses the package | Calls `h.traits` without `cor`/`cor.numr` | Closed-form two-sided test on the adjusted z statistic | Calls `h.traits` without `cor`/`cor.numr` | General one-trait calculation | Published closed form; `SINGLE_TRAIT_AFTER_SCREEN` |
| More than 16 traits in one direction | Uses the package | Stops the function | Records code 6, skips ASSET and continues the SNP loop | Prints a message, then continues toward ASSET | No fastASSET guard | Records, excludes and continues as in the paper script |
| Direction boundary at zero | Uses the package guard | Guard uses `>0` / `<=0` | Guard uses `>0` / `<=0` | Guard uses `>0` / `<=0` | Actual split uses `>=0` / `<0` | Uses `>=0` / `<0`; difference is confined to exact zero |
| Screening threshold | Example uses 0.05 | Function argument is accepted but 0.05 is passed internally | Fixed at 0.05 | Function argument is accepted but 0.05 is passed internally | No fastASSET screening | Honors the CLI value; default remains 0.05 |
| Meta | Not requested explicitly | Does not request the overall Meta component | Does not calculate the overall Meta component | Does not request the overall Meta component | `meta=TRUE` calculates Meta separately from subset searches | Optional, default TRUE; does not alter subset results |

The paper helper `ASSET_block_orth.R` is not used as production code. Its block
screening function captures `Neff` from the global environment instead of a
formal argument, and its use of `sapply` can simplify equal-length block lists
to a matrix. The standalone package corrected the first interface problem by
making `Neff` explicit and uses list-preserving block construction.

Two additional safety differences are intentional. fastASSETcli symmetrizes
the LDSC intercept matrix before normalization, and rejects a normalized matrix
whose minimum eigenvalue is not greater than the configurable tolerance
(`1e-8` by default). The paper script relies on downstream matrix operations to
fail instead of validating these conditions explicitly.

## Primary settings for this analysis

- Traits at input: 200–250.
- Input: one manifest containing one per-trait TSV/VCF file; no separate wide
  input is accepted.
- Quantitative NEF: total analyzed sample size, used directly.
- CPUs: `--ncores 90` on Linux/macOS.
- Primary screening threshold: `--scr-pthr 0.05`.
- Primary direction guard: `--max-traits-per-side 16`.
- A SNP exceeding that guard is written to
  `*_direction_limit_exceeded.tsv.gz`, skipped for ASSET enumeration, and does
  not stop analysis of other SNPs (unless the user explicitly enables
  `--fail-on-high TRUE` for the final exit policy).
- Meta output: `--include-meta TRUE` (also the default).

## Reference behavior used

- Upstream fastASSET standardizes as `beta / SE / sqrt(Neff)`, uses
  `SE = 1 / sqrt(Neff)`, pre-screens, and calls `h.traits` with `side=2` and
  `search=2`.
- ASSET `h.traits2` searches the non-negative and negative adjusted-effect
  groups separately. The wrapper guards in all three upstream code artifacts
  instead use `>0` / `<=0`; fastASSETcli follows the split actually enumerated
  by `h.traits2`.
- ASSET `meta=TRUE` returns `Meta = list(pval, beta, sd)`.
- GenomicSEM `ldsc()` returns a list containing `I` and `S`; trait names are
  assigned to `S`, so `S` is a necessary fallback when `I` has no usable names.

## Sources

- fastASSET source and tutorial: <https://github.com/gqi/fastASSET>
- Released 116-trait analysis: <https://github.com/gqi/fastASSET/blob/main/analysis/01_ASSET_GRASP_analysis.R>
- Paper block-screening helper: <https://github.com/gqi/fastASSET/blob/main/support_files/ASSET_block_orth.R>
- ASSET source: <https://github.com/sbstatgen/ASSET>
- ASSET Bioconductor manual: <https://bioconductor.org/packages/release/bioc/manuals/ASSET/man/ASSET.pdf>
- GenomicSEM LDSC source: <https://github.com/GenomicSEM/GenomicSEM/blob/master/R/ldsc.R>
- GenomicSEM munge source: <https://github.com/GenomicSEM/GenomicSEM/blob/master/R/munge_main.R>
- bcftools query documentation: <https://samtools.github.io/bcftools/bcftools.html#query>
- GWAS-VCF specification: <https://github.com/MRCIEU/gwas-vcf-specification>
- Qi et al. 2024: <https://www.nature.com/articles/s41467-024-51075-5>
