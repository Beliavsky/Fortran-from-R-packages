# kofnGA-fortran

Modern Fortran/FPM translation of the computational core of the R package
`kofnGA` 1.3 by Mark A. Wolters.

The library searches for a subset of exactly `k` indices from `1:n` that
minimizes a user-supplied objective function. It preserves the upstream GA
structure:

- rank-weighted tournament selection;
- crossover by sampling `k` elements from the union of the two parents;
- fixed-cardinality mutation by replacing selected genes with unused indices;
- optional elitism;
- reusable initial populations;
- best/minimum and population-average histories; and
- final population sorting and summary statistics.

## Build

```text
fpm build
fpm test
fpm run --example smallest_subset
```

The package is entirely free-form Fortran and uses:

```toml
[fortran]
implicit-typing = false
implicit-external = false
source-form = "free"
```

No BLAS/LAPACK or other external numerical library is required.

## Basic use

The objective is a normal Fortran procedure accepting an integer subset:

```fortran
function objective(subset) result(value)
  use kofnga, only : dp
  integer, intent(in) :: subset(:)
  real(dp) :: value
  value = sum(cost(subset))
end function objective
```

Then call:

```fortran
call kofn_ga(n, k, objective, result, control)
```

`kofnga_result` contains `bestsol`, `bestobj`, the final population/objective
values, and the complete best/average objective histories.

## Deliberate interface differences from R

The R-only `parallel`/`bigmemory` machinery is not translated. It changes how
objective evaluations are dispatched, not the genetic algorithm itself.
Fortran callers can parallelize an expensive objective internally if desired.

The Fortran implementation uses an explicitly specified Park-Miller RNG so a
seed is reproducible across Fortran builds. Therefore a Fortran seed is not
expected to reproduce R's `sample`, `rbinom`, and random-tie stream bit for bit.
The statistical selection/crossover/mutation rules are preserved.

See `API_MAPPING.md`, `ALGORITHM_NOTES.md`, and `VALIDATION.md`.
