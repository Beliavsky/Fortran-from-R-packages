#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript was not found on PATH." >&2
    exit 1
fi

if ! command -v fpm >/dev/null 2>&1; then
    echo "ERROR: fpm was not found on PATH." >&2
    exit 1
fi

r_results="$script_dir/r_results.csv"
fortran_results="$script_dir/fortran_results.csv"

echo "Fitting the selected GARCH models with R rugarch..."
Rscript reference.R "$r_results"

echo
echo "Fitting the selected GARCH models with the Fortran translation..."
fpm run --profile release -- "$fortran_results"

echo
echo "Both implementations completed successfully."
echo "R results:        $r_results"
echo "Fortran results:  $fortran_results"
