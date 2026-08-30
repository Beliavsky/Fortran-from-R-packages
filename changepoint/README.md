# changepoint — modern Fortran computational translation

This package translates the computational algorithms of R package
`changepoint` 2.3 to modern free-form Fortran. It targets the root of
`Beliavsky/Fortran-from-R-packages` and is designed to build with the Fortran
Package Manager (FPM).

The translation covers likelihood-based AMOC, pruned PELT, Binary
Segmentation, exact Segment Neighborhood, CUSUM/CSS nonparametric methods,
CROPS, Normal regression changepoints, penalty calculations, asymptotic AMOC
transforms, and segment parameter estimation. R-specific S4, plotting,
formula, and time-series presentation code is deliberately omitted. See
`API_COVERAGE.md` for the detailed mapping.

## Dependencies

Place this directory beside the shared packages in the repository root:

```text
Fortran-from-R-packages/
├── changepoint/
├── rfortran-core/
└── rfortran-linalg/
```

`fpm.toml` uses sibling path dependencies. No dependency sources are vendored.
`dp` is imported from `rfortran-core`'s `r_kinds` module and re-exported by the
public `changepoint` module. All maintained Fortran code uses `real(dp)` and
`_dp` real literals.

## Build and test

```console
cd changepoint
fpm build
fpm test
fpm run --example mean_changes
fpm run --example regression_change
```

For a release-tree audit, run:

```console
python tools/release_check.py
```

The release checker performs FPM build/test, source-hygiene checks, and
`fpm clean --all`.

## Basic PELT example

```fortran
program demo
use changepoint, only : dp, changepoint_result
use changepoint, only : cp_pelt, cp_cost_mean_normal
implicit none
real(dp) :: x(12)
type(changepoint_result) :: result

x = [0.0_dp, 0.1_dp, -0.1_dp, 0.05_dp, &
     5.0_dp, 5.1_dp, 4.9_dp, 5.05_dp, &
    -2.0_dp, -2.1_dp, -1.9_dp, -2.05_dp]

call cp_pelt(x, cp_cost_mean_normal, 1.0_dp, 2, result)
if (result%status == 0) then
    print *, result%cpts
end if
end program demo
```

This prints changepoints `4 8`. As in upstream `changepoint`, a changepoint is
the **last observation of the preceding segment**.

## Cost models

The public cost codes are:

- `cp_cost_mean_normal`
- `cp_cost_var_normal`
- `cp_cost_meanvar_normal`
- `cp_cost_exponential`
- `cp_cost_gamma`
- `cp_cost_poisson`

Gamma accepts a `shape` argument. Normal-variance methods can accept a
`known_mean`. Exponential observations must be nonnegative; Gamma observations
must be strictly positive; Poisson observations must be nonnegative integers.

## Regression changepoints

`cp_regression_amoc` and `cp_regression_pelt` take the response vector and the
design matrix directly. Intercepts are therefore explicit columns of the design
matrix rather than being inserted by R formula machinery. Least squares uses
`rfortran-linalg`'s SVD solver.

## Licensing and provenance

The upstream package declares `License: GPL` without a version qualifier.
Original package metadata and citation information are retained under
`upstream/`; see `NOTICE.md`, `UPSTREAM.md`, and `LICENSE.note`.
