# isotone-fortran

Modern free-form Fortran/FPM translation of the computational code in the R
package `isotone` 1.1-2.

The package implements generalized pool-adjacent-violators regression (PAVA),
a primal active-set method for arbitrary pairwise order restrictions, the
predefined isotone loss functions, and regression with linear inequality
restrictions on fitted values.

## Main API

```fortran
use isotone

type(gpava_result) :: pfit
type(active_set_result) :: afit
type(active_set_options) :: control

call gpava_fit(y, pfit, weights=w)

control%solver = ISO_LS
allocate(control%weights(size(y)))
control%weights = w
call active_set(order_pairs, y, afit, control)
```

`order_pairs(k,:) = [i,j]` means fitted value `j` must be greater than or
equal to fitted value `i`, matching the R package convention.

The repeated-measurement PAVA API is `gpava_fit_repeated` with an `n x r`
response matrix.

Individual upstream solver functions are available as `ls_solver`,
`d_solver`, `p_solver`, `lf_solver`, `s_solver`, `o_solver`, `a_solver`,
`e_solver`, `h_solver`, `i_solver`, `m_solver`, and callback-driven
`f_solver`.

The NNLS-based functions are:

```fortran
call mregnn(x, y, a, result)
call mregnn_monotone(x, y, result)
call mregnn_positive(x, y, result)
```

## Implemented functionality

- generalized PAVA with weighted mean, weighted median, and weighted fractile
  block solvers;
- primary, secondary, and tertiary tie handling;
- repeated-measurement PAVA;
- arbitrary pairwise partial orders in the primal active-set algorithm;
- weighted least squares and full generalized least squares;
- L1, quantile, Chebyshev, Poisson, Lp, asymmetric LS, smoothed L1, Huber,
  and SILF losses;
- user-defined differentiable convex losses through a Fortran callback;
- KKT multipliers and stationarity/primal/dual/complementary-slackness
  diagnostics;
- `mregnn`, monotone `mregnnM`, and positive `mregnnP` equivalents using the
  supplied `nnls-fortran` translation.

Plotting, S3 printing/summary methods, R formula/list dispatch, and bundled
datasets are not part of the numerical library.

## Build

```text
fpm build
fpm test
```

The manifest uses only free source form and disables implicit typing and
implicit external interfaces. No BLAS/LAPACK dependency is required.

## Licensing

The upstream `isotone` package is GPL-2. The supplied NNLS Fortran dependency
is GPL-2.0-or-later and is vendored in `src/vendor_nnls`. The combined package
is distributed under GPL-2. See `COPYING`, `COPYRIGHTS`, `NNLS_NOTICE.md`, and
`UPSTREAM.md`.
