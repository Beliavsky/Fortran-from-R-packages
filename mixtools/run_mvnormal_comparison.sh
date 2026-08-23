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

data_file="example/mvnormal_mixture_data.txt"

if [[ ! -f "$data_file" ]]; then
    echo "Generating multivariate-normal mixture data..."
    Rscript "example/generate_mvnormal_mixture.R"
else
    echo "Using existing data file \"$data_file\"."
fi

echo
echo "Fitting the mixture with R mixtools..."
Rscript "example/fit_mvnormal_mixture.R"

echo
echo "Fitting the mixture with the Fortran translation..."
fpm run --example fit_mvnormal_mixture

echo
echo "Comparison completed successfully."
echo "Reports:"
echo "  example/mvnormal_mixture_fit_r.txt"
echo "  example/mvnormal_mixture_fit_fortran.txt"
