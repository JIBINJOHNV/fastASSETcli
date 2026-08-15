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
