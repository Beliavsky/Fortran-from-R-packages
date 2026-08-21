# Porting notes

## Architecture

The R package is translated into small numerical modules rather than emulating R
objects. Array arguments replace R vectors/data frames, and derived types replace
list/data-frame results. The project has no runtime dependency on R, MASS, or boot.

The only use of MASS in the upstream computational source is plotting (`eqscplot`).
The only use of boot is the ordinary bootstrap inside `vm.bootstrap.ci`; the Fortran
port implements that specialized bootstrap directly.

## Modified Bessel functions

`dvm`, `A1`, and related routines use scaled modified-Bessel approximations so that
large concentration values do not overflow through `exp(kappa)`. `pvm` is evaluated
by adaptive Simpson integration of the same von Mises density rather than by the
upstream infinite Fourier series. This is mathematically equivalent and avoids
unscaled high-order Bessel overflow.

## Preserved Rao table

`data/rao.table.rda` contains a 43 x 4 matrix of critical values. The matrix was
decoded from the upstream R serialization and embedded verbatim in
`src/circstats_rao_table.f90`.

## Upstream fixes / numerical regularizations

1. **`rstable`, skewed alpha=1 branch.** The R source references an undefined symbol
   `c` in `log(c)`. The surrounding GSL-derived formula and the function's `scale`
   argument make the intended term the scale adjustment. The port therefore uses
   `log(scale)`.

2. **`rcard` with negative concentration.** Documentation permits `abs(r) < 0.5`,
   but the upstream rejection envelope `(1 + 2*r)/(2*pi)` is not an upper bound when
   `r < 0`. The port uses `(1 + 2*abs(r))/(2*pi)`.

3. **Uniform limiting cases.** `rvm(kappa=0)`, `rtri(r=0)`, and wrapped-normal/Cauchy
   zero-concentration limits are handled explicitly rather than allowing divisions
   by zero in the literal formulas.

4. **`pvm(2*pi)`.** The port returns 1 at the full-circle endpoint, which is the
   cumulative-probability limit. The R implementation first applies modulo and thus
   maps exactly `2*pi` back to zero.

5. **Iteration safeguards.** Wrapped-Cauchy MLE and density summations have explicit
   iteration limits so invalid inputs cannot create unbounded loops.

## Circular regression

The upstream regression constructs dense n x n projection matrices. The Fortran
implementation uses algebraically equivalent residualized regressors for the
higher-order test, reducing memory while producing the same statistics.

## Graphics

No rendering code is included. `pp_fit` retains the non-graphics return value from
`pp.plot` (von Mises mean direction and concentration estimate).
