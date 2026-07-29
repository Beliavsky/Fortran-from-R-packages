# mcrp-fortran

A modern Fortran 2018 translation of the computational algorithms in the R
package **mcrp 0.0-1**, "Multiple criteria risk parity optimization".

The library computes centered second-, third-, and fourth-order co-moments,
portfolio variance/skewness/kurtosis decompositions, and portfolios that balance
risk contributions across any selected combination of those three criteria.
It is self-contained and has no external numerical-library dependency.

## Implemented algorithms

- `M2`, `M3`, and `M4` centered co-moment matrices/tensors
- `pm2`, `pm3`, and `pm4` portfolio centered moments
- `dm2`, `dm3`, and `dm4` moment derivatives
- `cm2`, `cm3`, and `cm4` raw-moment contribution decompositions
- Portfolio variance, standardized skewness, and standardized kurtosis
- The original portfolio contribution decompositions
- Multiple-criteria risk-parity objective
- Bounded derivative-free optimization
- Optional long-only, box-bounded, or unrestricted raw parameters
- Typed optimization results and explicit status codes

All 22 exported numerical entry points from the R package are represented.

## Build with FPM

```sh
fpm test
fpm run --target mcrp_demo
fpm run --example moment_example
```

FPM was not installed in the translation environment, so the included manifest
was parsed independently and all targets were also built with GNU Fortran using:

```sh
./build_gfortran.sh strict
./build_gfortran.sh release
```

## Basic use

```fortran
use mcrp_module

real(dp) :: returns(100, 4)
real(dp) :: start(4), lower(4)
type(mcrp_result) :: fit

start = 0.25_dp
lower = 0.0_dp
call mcrp(start, returns, fit, lower=lower)

if (fit%converged) then
   print *, fit%weights
end if
```

To optimize only variance risk parity:

```fortran
call mcrp(start, returns, fit, lower=lower, &
   active=[.true., .false., .false.])
```

The final weights are normalized by their L1 norm, matching the R package:
`weights = raw_parameters / sum(abs(raw_parameters))`.

## Important compatibility details

The source package labels `PortSkewDeriv` and `PortKurtDeriv` as derivatives,
but their formulas are not the ordinary gradients of standardized skewness and
kurtosis. They are decomposition vectors constructed so that weighted
components sum to the reported portfolio measure. The Fortran port preserves
those formulas exactly.

The source optimizer is R's `nlminb`. This port uses a self-contained bounded
Nelder-Mead implementation. The objective and normalization are preserved, but
iteration paths and final local solutions can differ.

See `PORTING.md` for all compatibility details.

## License

GPL-3.0-only, matching the original package metadata. The original package tree
is retained under `original/` for provenance.
