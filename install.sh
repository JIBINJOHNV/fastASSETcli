#!/usr/bin/env bash
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
package_dir="$script_dir"
r_library="${FASTASSET_R_LIB:-${HOME}/.local/share/fastassetcli/R-library}"
bin_dir="${FASTASSET_BIN_DIR:-${HOME}/.local/bin}"

if ! command -v R >/dev/null 2>&1; then
  echo "ERROR: R is not installed or is not on PATH." >&2
  exit 2
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript is not installed or is not on PATH." >&2
  exit 2
fi

mkdir -p "$r_library" "$bin_dir"

if [[ "${FASTASSET_SKIP_DEPENDENCIES:-0}" != "1" ]]; then
  Rscript "${script_dir}/tools/install_dependencies.R" "$r_library"
fi

R_LIBS_USER="$r_library" R CMD INSTALL \
  --library="$r_library" \
  "$package_dir"

installed_launcher=$(R_LIBS_USER="$r_library" Rscript -e \
  'cat(system.file("exec", "fastasset", package = "fastASSETcli"))')
if [[ -z "$installed_launcher" || ! -f "$installed_launcher" ]]; then
  echo "ERROR: The installed fastasset launcher was not found." >&2
  exit 2
fi

chmod +x "$installed_launcher"
ln -sfn "$installed_launcher" "${bin_dir}/fastasset"

echo "Installed fastASSETcli into: $r_library"
echo "Installed shell command: ${bin_dir}/fastasset"
if [[ ":${PATH}:" != *":${bin_dir}:"* ]]; then
  echo "Add this directory to PATH before using the command:"
  echo "  export PATH=\"${bin_dir}:\$PATH\""
fi
echo "Test the installation with:"
echo "  fastasset --version"
echo "  fastasset --help"
