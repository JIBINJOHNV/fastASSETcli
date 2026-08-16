# fastASSETcli 0.4.9

- Name converted VCF/BCF tables and metadata with both the zero-padded manifest
  order and a filesystem-safe manifest trait label, for example
  `trait_0002_ANXA1_P04083_OID30617_v1_Inflammation_II.tsv.gz`.
- Apply the same informative stem to in-progress `bcftools-query` and stderr
  files, so every parallel worker can be linked to its trait while running.
- Updated the prepared-input cache signature so an earlier cache containing
  numeric-only VCF-conversion filenames is not reused.

# fastASSETcli 0.4.8

- Import `data.table` in `NAMESPACE`, making the installed package explicitly
  data.table-aware. This fixes the `cedta()` runtime error raised when package
  functions use `:=`, `.N`, or grouped `data.table` expressions.
- Added a regression assertion that the installed package namespace contains
  the required `data.table` import.

# fastASSETcli 0.4.7

- In existing-LDSC mode, compare manifest traits with the native GenomicSEM
  trait order from `colnames(LDSCoutput$S)` before processing summary files.
- Keep only common traits, preserve their manifest order, subset manifest rows
  and VCF SAMPLE values consistently, and reorder both `$I` axes accordingly.
  Excluded rows are removed before their file paths or samples are validated.
- Warn separately for traits missing from LDSC and traits missing from the
  manifest, and write `manifest_ldsc_trait_match.tsv` with the complete audit.
- Write every `DIRECTION_LIMIT_EXCEEDED` SNP to the dedicated
  `*_direction_limit_exceeded.tsv.gz` file while continuing with other SNPs.
- Added regression tests for common-trait ordering, sample-row subsetting and
  dedicated direction-limit records.
- No additional scientific QC filters were added; upstream-QC-completed input
  is used subject only to structural requirements needed for valid analysis.

# fastASSETcli 0.4.6

- Made every temporary `bcftools query` table self-describing by writing one
  canonical tab-separated header before appending bcftools's data stream.
  Interrupted conversions can now be inspected with ordinary tools such as
  `head` without relying on positional knowledge.
- Added a dedicated header-aware query reader that validates both the field
  count and canonical names. The large data stream is written only once; no
  second pass or full-file rewrite is required.

# fastASSETcli 0.4.5

- Fixed quantitative and binary VCF conversion when `data.table::fread()` is
  configured by the R environment to return a base `data.frame`. Valid-row
  filtering is now explicitly row-oriented and no longer fails with
  `undefined columns selected`.
- Made the bcftools query reader explicitly return a `data.table`, validate the
  extracted field count before assigning names, and report the expected VCF
  fields in any future column-count error.

# fastASSETcli 0.4.4

- Reorganized the main README into a user-first sequence covering
  installation, manifest creation, LDSC mode selection, scientific input
  rules, recommended settings and primary outputs.
- Added a linked complete workflow guide with a full pipeline flowchart,
  fifteen processing steps, all decision routes, QC statuses, resume behavior
  and output-directory structure.
- Added clearer copy-ready commands for existing-LDSC and automatic-LDSC
  modes and consolidated the quantitative, binary, tabular and VCF sample-size
  rules into comparison tables.

# fastASSETcli 0.4.3

- Changed incompatible cross-trait allele handling from a fatal error to an
  auditable SNP-trait exclusion. BETA, SE, and NEF are set to `NA` only for
  the incompatible trait at that SNP; compatible traits at the SNP remain.
- Added `fastasset_failed_allele_alignments.tsv` with the trait, SNP ID,
  source file, input row, reference alleles, incoming alleles, reason, and
  action for every excluded SNP-trait observation.
- Confirmed through regression tests that incompatible observations are not
  passed to ASSET. The existing per-SNP validity filter analyzes the remaining
  traits, or records `INSUFFICIENT_VALID_TRAITS` when fewer than the configured
  minimum remain.
- Updated the input-preparation cache signature and completion requirements so
  a prepared input is reusable only when the failed-alignment report exists.

# fastASSETcli 0.4.2

- Retained the efficient single SNP-ID index plus exact/swapped allele
  comparison used to construct the FastASSET wide table; no repeated physical
  merges are performed.
- Expanded `fastasset_input_build_audit.tsv` with exact matches, swapped
  matches, reference-establishing SNPs, final union size, SNPs absent from each
  trait, and incompatible-match counts. Incompatible pairs remain fatal.
- Added regression coverage for exact matches, sign-corrected swapped matches,
  missing union SNPs, and rejection of strand-complement/incompatible pairs.
- Updated the input-preparation cache signature for the expanded audit.

# fastASSETcli 0.4.1

- Binary VCF LDSC continues to use per-SNP total sample size
  `N=NC+NCO`, paired with one trait-level sample prevalence calculated as
  `median(NC)/(median(NC)+median(NCO))` across valid SNP rows.
- Binary FastASSET continues to use the original per-SNP definition
  `NC*NCO/(NC+NCO)`. VCF `NEF` is not substituted for either binary rule.
- VCF header `TotalCases` and `TotalControls` are retained as provenance but
  no longer determine binary LDSC sample prevalence.
- Updated conversion provenance and cache signatures for the prevalence rule.

# fastASSETcli 0.4.0

- Removed the public `--fastasset-input` option. The per-trait
  `--sumstats-manifest` is now the only summary-statistic input in both
  existing-LDSC and automatic-LDSC modes.
- Added deterministic construction of the internal `ID` plus
  `<trait>.Beta/.SE/.NEF` table in manifest row order.
- Added cross-trait A1/A2 verification: swapped pairs flip BETA and
  incompatible pairs fail before analysis.
- Added original-reference binary FastASSET effective sample size,
  `Ncase*Ncontrol/(Ncase+Ncontrol)`, while preserving direct quantitative NEF.
- Added preparation provenance, allele-reference and per-trait build-audit
  outputs with signature-protected reuse.

# fastASSETcli 0.3.0

- Added automatic `bcftools` conversion for `.vcf`, `.vcf.gz`, `.vcf.bgz`,
  and `.bcf` files listed in the LDSC manifest.
- Added explicit GWAS-VCF mapping: `A1=ALT`, `A2=REF`, `BETA=ES`, `INFO=SI`,
  `MAF=min(AF,1-AF)`, and `P=10^(-LP)` with positive-double protection.
- Added binary-trait handling when `POPULATION_PREV` is supplied: per-SNP
  `N=NC+NCO` and automatic `SAMPLE_PREV=NC/(NC+NCO)` derivation.
- Added quantitative VCF handling when both prevalences are absent: per-SNP
  `FORMAT/NEF` is used directly as total sample size.
- Added VCF header/FORMAT validation, multi-sample selection, resumable
  conversion outputs, conversion metadata, `--bcftools`, and `--vcf-cores`.

# fastASSETcli 0.2.0

- Made `--ldsc-rdata` optional by adding automatic GenomicSEM munging and LDSC.
- Added a validated per-trait summary-statistic manifest with exact FastASSET
  trait-set matching and deterministic reordering.
- Added CLI options for HapMap3 and LD-score references, prevalence, optional
  sample-size overrides, column mapping, munging workers, chromosomes,
  jackknife blocks, and chi-square filtering.
- Added separate, signature-addressed LDSC output directories containing the
  generated `.RData`, munged files, logs, a resolved manifest, and provenance.
- Preserved the existing audited `.RData`/`.rda`/`.rds` loading mode.

# fastASSETcli 0.1.0

- Added an installable R package and `fastasset` shell command.
- Added CLI-only execution with no configuration file.
- Added separate GenomicSEM `.RData`, `.rda`, and `.rds` loading.
- Added audited trait-name recovery and exact two-axis LDSC matrix reordering.
- Added LDSC covariance-to-correlation normalization and matrix validation.
- Added reference-aligned `side=2`, `search=2` fastASSET execution.
- Added ASSET fixed-effect Meta output for screened, selection-adjusted traits.
- Added per-SNP status, severity, directional counts, and resumable chunks.
