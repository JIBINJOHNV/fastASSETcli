# Installable fastASSET CLI package

This repository installs the reviewed pipeline as the R package `fastASSETcli` and
provides the shell command `fastasset`. It does not install under the name
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
2. installs `data.table`, Bioconductor `ASSET`, and their dependencies when
   needed;
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

Set `FASTASSET_SKIP_DEPENDENCIES=1` only when `ASSET` and `data.table` are
already visible to the selected R installation.

## Run the analysis

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

## Manual R package installation

If dependencies are already installed, the package directory can also be
installed normally:

```bash
R CMD INSTALL .
# or
R CMD INSTALL release/fastASSETcli_0.1.0.tar.gz
```

The executable inside an installed package can be located with:

```bash
Rscript -e 'cat(system.file("exec", "fastasset", package="fastASSETcli"))'
```

## Scientific behavior

- Quantitative `.NEF` is used directly as total analyzed sample size.
- GenomicSEM `$I` is loaded from a separate `.RData`, `.rda`, or `.rds` file.
- Trait names are recovered from `I`, `S`, or `S_Stand` and reordered to the
  `<trait>.Beta` column order.
- The LDSC covariance is normalized as `D^(-1/2) I D^(-1/2)`.
- ASSET is called with `side=2`, `search=2`, and `meta=TRUE` by default.
- The primary guard remains 16 screened traits per direction.
- Meta output is across screened, selection-adjusted traits.

The detailed supplied-code comparison is installed with the package under
`docs/REFERENCE_COMPARISON.md`.

## Release archive

`release/fastASSETcli_0.1.0.tar.gz` is the installable source-package archive
for version 0.1.0. The editable package source remains at the repository root;
the archive is not a substitute for version-controlled source code.

## Upstream terms

This workflow contains adapted fastASSET logic and calls ASSET. Review the
terms of both upstream projects before public redistribution. The upstream
fastASSET repository does not currently declare a standard license in its
`DESCRIPTION`, while ASSET provides separate NCI terms. Keeping a new GitHub
repository private is recommended until redistribution permission is clear.
