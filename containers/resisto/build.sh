#!/usr/bin/env bash
# Build the resisto-GEM Apptainer image with the CPLEX Python API from a licensed CPLEX Studio.
#
#   bash containers/resisto/build.sh /path/to/CPLEX_Studio2211 [output.sif]
#
# The CPLEX Python API must match the image's Python (3.10). The script looks for the directory
# holding CPLEX's setup.py for Python 3.10 on x86-64 Linux inside the CPLEX Studio folder; if your
# release lays this out differently, copy that directory to containers/resisto/cplex_python yourself.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
studio="${1:?usage: build.sh /path/to/CPLEX_Studio [output.sif]}"
out="${2:-resisto-gem.sif}"

candidates=$(find "$studio" -path '*python*' -name setup.py 2>/dev/null | grep -E '3\.10|/python/setup.py$' || true)
if [[ -z "$candidates" ]]; then
    echo "No CPLEX Python setup.py for Python 3.10 found under $studio" >&2
    exit 1
fi
src=$(dirname "$(echo "$candidates" | grep -E 'x86-64_linux' | head -n 1 || echo "$candidates" | head -n 1)")
echo "Using CPLEX Python API from: $src"

rm -rf "$here/cplex_python"
cp -r "$src" "$here/cplex_python"
cd "$here"
apptainer build "$out" resisto.def
echo "Built $out. Keep it private: it contains licensed CPLEX files."
