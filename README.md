# Installable fastASSET CLI package

This repository installs the reviewed pipeline as the R package `fastASSETcli` and
provides the shell command `fastasset`. Version 0.4.0 uses one per-trait
summary-statistic manifest as the complete analysis input. It can either load an
existing GenomicSEM LDSC object or generate the object automatically from
tabular or GWAS-VCF summary statistics before running FastASSET. The multi-trait
wide table required internally by FastASSET is constructed automatically; there
is no `--fastasset-input` option. The package does not install under the name
`fastASSET`, so it does not overwrite or shadow the upstream project.

## One-command installation

From this directory:

```bash
chmod +x install.sh
./install.sh
```

The installer:

1. creates an isolated user R library under
   `~/.local/share/fastassetcli/R-library`;
2. installs `data.table`, Bioconductor `ASSET`, GenomicSEM, and their
   dependencies when needed;
3. installs the `fastASSETcli` R package; and
4. links the executable to `~/.local/bin/fastasset`.

If `~/.local/bin` is not already on `PATH`, add it:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Then verify:

```bash
fastasset --version
fastasset --help
```

The installation locations can be changed without a configuration file:

```bash
FASTASSET_R_LIB=/opt/R/fastasset \
FASTASSET_BIN_DIR=/opt/bin \
./install.sh
```

Set `FASTASSET_SKIP_DEPENDENCIES=1` only when `ASSET`, `data.table`, and
`GenomicSEM` are already visible to the selected R installation.

## Run with an existing LDSC object

```bash
fastasset \
  --sumstats-manifest /data/manifest.tsv \
  --ldsc-rdata /data/all_traits_LDSCoutput.RData \
  --ldsc-object-name LDSCoutput \
  --output-dir /results/fastasset \
  --run-name quantitative_250_traits \
  --ncores 90 \
  --chunk-size 1000 \
  --scr-pthr 0.05 \
  --max-traits-per-side 16 \
  --include-meta TRUE
```

Everything is passed through the shell CLI; no YAML or other configuration
file is read.

## Generate LDSC automatically

When `--ldsc-rdata` is omitted, the CLI runs `GenomicSEM::munge()`, runs
`GenomicSEM::ldsc()`, saves the complete GenomicSEM object separately, and then
continues with the same audited FastASSET workflow.

Create a tab- or comma-delimited manifest with exactly one row per trait. Each
`FILE` may be an ordinary summary-statistic table or a `.vcf`, `.vcf.gz`,
`.vcf.bgz`, or `.bcf` file:

```text
TRAIT    FILE                              N    SAMPLE    SAMPLE_PREV    POPULATION_PREV
trait_1  /data/trait_1_sumstats.tsv.gz     NA   NA        NA             NA
trait_2  /data/trait_2.vcf.gz              NA   NA        NA             NA
```

The manifest row order is the canonical trait order for both the internally
constructed FastASSET table and LDSC. When an existing LDSC object is supplied,
both axes of its intercept matrix are reordered to this manifest order.
`NAME`, `SPREV`, and `PPREV` are accepted as aliases for compatibility with the
supplied pipeline.

For quantitative traits, set both prevalence fields to `NA`. `N` is optional:
leave it `NA` when each raw summary-statistic file already contains a valid
sample-size column. A manifest `N` is a constant per-trait sample size used by
both internal FastASSET input construction and GenomicSEM munging when a
per-SNP sample-size column is absent.

For VCF/BCF inputs, leave manifest `N` as `NA`. A one-sample VCF is selected
automatically. If a VCF has more than one sample, set `SAMPLE` to the exact VCF
sample ID. VCF conversion applies these audited rules:

| Input situation | LDSC sample size | Prevalence handling |
|---|---|---|
| `POPULATION_PREV` supplied | Per-SNP `N = FORMAT/NC + FORMAT/NCO` | Binary trait; `SAMPLE_PREV` is calculated as `NC/(NC+NCO)` using the unambiguous maximum-total-sample row. Header `TotalCases`/`TotalControls` takes priority when present. |
| `POPULATION_PREV` is `NA` | Per-SNP `N = FORMAT/NEF` | Quantitative trait; both prevalences remain `NA`. Here `NEF` is used directly as total analyzed sample size. |

The supplied GWAS-VCF convention is handled explicitly: `BETA=FORMAT/ES`,
`A1=ALT`, `A2=REF`, `MAF=min(FORMAT/AF, 1-FORMAT/AF)`,
`INFO=FORMAT/SI` when present, and raw `P=10^(-FORMAT/LP)`. Extremely small P values that
cannot be represented by R are capped at `.Machine$double.xmin`, never changed
to zero. A VCF declaring `StudyType=CaseControl` without
`POPULATION_PREV` is rejected.

For every TSV trait file, input construction identifies SNP ID, A1, A2, BETA,
SE and sample size. Quantitative `NEF` is used directly as total analyzed N.
For binary traits, the original fastASSET definition is used:
`Neff = NCASE*NCONTROL/(NCASE+NCONTROL)`. If case/control columns are absent,
the equivalent `N*SAMPLE_PREV*(1-SAMPLE_PREV)` is used. Effect alleles are
checked across all traits. An exactly swapped allele pair causes BETA to be
flipped; an incompatible pair is a fatal harmonization error.

```bash
fastasset \
  --sumstats-manifest /data/manifest.tsv \
  --hm3 /refs/w_hm3.snplist \
  --ld-ref /refs/eur_w_ld_chr \
  --output-dir /results/fastasset \
  --run-name quantitative_250_traits \
  --ncores 90 \
  --munge-cores 90 \
  --vcf-cores 90 \
  --chunk-size 1000 \
  --scr-pthr 0.05 \
  --max-traits-per-side 16 \
  --include-meta TRUE
```

If the input columns use nonstandard names, pass a GenomicSEM mapping entirely
through the CLI. For example:

```bash
--munge-column-map 'SNP=ID,A1=EA,A2=NEA,effect=BETA,SE=SE,N=NEF,P=P'
```

Use `--wld-ref` only when the LD-score weights are stored separately from
`--ld-ref`. By default, the same directory is used for both.

`bcftools` must be on `PATH` for VCF inputs. Alternatively provide it through
the CLI, for example `--bcftools /opt/bcftools/bin/bcftools`. This dependency
is not required when every manifest `FILE` is already tabular.

Automatic outputs are written under:

```text
<output-dir>/<run-name>_prepared_input_<signature>/
├── <run-name>_fastasset_wide.tsv.gz
├── fastasset_allele_reference.tsv.gz
├── fastasset_input_build_audit.tsv
├── manifest_resolved.tsv
├── input_preparation_provenance.tsv
└── vcf_converted/                    # present only for VCF/BCF inputs

<output-dir>/<run-name>_generated_ldsc_<signature>/
├── <run-name>_LDSCoutput.RData
├── ldsc_manifest_resolved.tsv
├── ldsc_generation_provenance.tsv
├── <run-name>_ldsc.log
└── munged/
```

The signature includes the trait order, summary-statistic metadata, reference
metadata, prevalence values, and LDSC/munging settings. A completed matching
object is reused; incompatible or incomplete outputs are not silently reused.

Important for 200–250 traits: current GenomicSEM automatically increases the
jackknife blocks above 18 traits to
`((number_of_traits + 1) * (number_of_traits + 2) / 2) + 1`. Munging uses
`--munge-cores`, but `GenomicSEM::ldsc()` itself does not expose a worker-count
argument. The CLI reports this before LDSC starts.

## Manual R package installation

If dependencies are already installed, the package directory can also be
installed normally:

```bash
R CMD INSTALL .
# or
R CMD INSTALL release/fastASSETcli_0.4.0.tar.gz
```

The executable inside an installed package can be located with:

```bash
Rscript -e 'cat(system.file("exec", "fastasset", package="fastASSETcli"))'
```

## Scientific behavior

- The manifest is the only summary-statistic input; its row order defines the
  trait order, and the FastASSET wide table is generated internally.
- Quantitative `NEF` is used directly as total analyzed sample size.
- Binary fastASSET `Neff` is calculated as
  `NCASE*NCONTROL/(NCASE+NCONTROL)`, matching the upstream reference.
- Automatic LDSC accepts tabular or GWAS-VCF manifest files; VCF binary
  prevalence and sample size are derived from `NC/NCO`, while quantitative VCF
  `NEF` is used directly.
- GenomicSEM `$I` is loaded from an existing object or generated and saved as a
  separate `.RData` file.
- Trait names are recovered from `I`, `S`, or `S_Stand` and reordered to the
  manifest row order.
- The LDSC covariance is normalized as `D^(-1/2) I D^(-1/2)`.
- ASSET is called with `side=2`, `search=2`, and `meta=TRUE` by default.
- The primary guard remains 16 screened traits per direction.
- Meta output is across screened, selection-adjusted traits.

The detailed supplied-code comparison is installed with the package under
`docs/REFERENCE_COMPARISON.md`.

## Release archive

`release/fastASSETcli_0.4.0.tar.gz` is the current installable source-package
archive. Earlier archives remain available for reproducibility. The editable
package source remains at the repository root; an archive is not a substitute
for version-controlled source code.

## Upstream terms

This workflow contains adapted fastASSET logic and calls ASSET. Review the
terms of both upstream projects before public redistribution. The upstream
fastASSET repository does not currently declare a standard license in its
`DESCRIPTION`, while ASSET provides separate NCI terms. Keeping a new GitHub
repository private is recommended until redistribution permission is clear.
