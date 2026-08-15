# Installable fastASSET CLI package

This repository installs the reviewed pipeline as the R package `fastASSETcli` and
provides the shell command `fastasset`. Version 0.2.0 can either load an
existing GenomicSEM LDSC object or generate the object automatically before
running FastASSET. It does not install under the name
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
  --fastasset-input /data/fastasset_wide.tsv.gz \
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

Create a tab- or comma-delimited manifest with exactly one row for every
`<trait>.Beta` column in the FastASSET input:

```text
TRAIT    FILE                              N    SAMPLE_PREV    POPULATION_PREV
trait_1  /data/trait_1_sumstats.tsv.gz     NA   NA             NA
trait_2  /data/trait_2_sumstats.tsv.gz     NA   NA             NA
```

The trait names must match exactly. Manifest rows can be in any order because
the CLI validates and reorders them to the FastASSET `.Beta` column order.
`NAME`, `SPREV`, and `PPREV` are accepted as aliases for compatibility with the
supplied pipeline.

For quantitative traits, set both prevalence fields to `NA`. `N` is optional:
leave it `NA` when each raw summary-statistic file already contains a valid
sample-size column. A manifest `N` is a constant per-trait override used only
by GenomicSEM munging; it does not replace the per-SNP `<trait>.NEF` used by
FastASSET.

```bash
fastasset \
  --fastasset-input /data/fastasset_wide.tsv.gz \
  --sumstats-manifest /data/ldsc_manifest.tsv \
  --hm3 /refs/w_hm3.snplist \
  --ld-ref /refs/eur_w_ld_chr \
  --output-dir /results/fastasset \
  --run-name quantitative_250_traits \
  --ncores 90 \
  --munge-cores 90 \
  --chunk-size 1000 \
  --scr-pthr 0.05 \
  --max-traits-per-side 16 \
  --include-meta TRUE
```

If the input columns use nonstandard names, pass a GenomicSEM mapping entirely
through the CLI. For example:

```bash
--munge-column-map 'SNP=ID,A1=EA,A2=NEA,effect=BETA,N=NEF,P=P'
```

Use `--wld-ref` only when the LD-score weights are stored separately from
`--ld-ref`. By default, the same directory is used for both.

Automatic outputs are written under:

```text
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
R CMD INSTALL release/fastASSETcli_0.2.0.tar.gz
```

The executable inside an installed package can be located with:

```bash
Rscript -e 'cat(system.file("exec", "fastasset", package="fastASSETcli"))'
```

## Scientific behavior

- Quantitative `.NEF` is used directly as total analyzed sample size.
- GenomicSEM `$I` is loaded from an existing object or generated and saved as a
  separate `.RData` file.
- Trait names are recovered from `I`, `S`, or `S_Stand` and reordered to the
  `<trait>.Beta` column order.
- The LDSC covariance is normalized as `D^(-1/2) I D^(-1/2)`.
- ASSET is called with `side=2`, `search=2`, and `meta=TRUE` by default.
- The primary guard remains 16 screened traits per direction.
- Meta output is across screened, selection-adjusted traits.

The detailed supplied-code comparison is installed with the package under
`docs/REFERENCE_COMPARISON.md`.

## Release archive

`release/fastASSETcli_0.2.0.tar.gz` is the current installable source-package
archive. Version 0.1.0 remains available for reproducibility. The editable
package source remains at the repository root; an archive is not a substitute
for version-controlled source code.

## Upstream terms

This workflow contains adapted fastASSET logic and calls ASSET. Review the
terms of both upstream projects before public redistribution. The upstream
fastASSET repository does not currently declare a standard license in its
`DESCRIPTION`, while ASSET provides separate NCI terms. Keeping a new GitHub
repository private is recommended until redistribution permission is clear.
