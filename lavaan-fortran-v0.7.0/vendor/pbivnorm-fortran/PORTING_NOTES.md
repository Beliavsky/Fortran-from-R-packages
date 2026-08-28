# Porting notes

- The numerical core is a free-form Fortran 2018 translation of the upstream
  Genz/Ge bivariate-normal routine.
- Unlike the original fixed-form file, all calls have explicit module
  interfaces; the FPM project therefore keeps both implicit typing and implicit
  external procedures disabled.
- `pbivnorm` is elemental, which provides natural scalar broadcasting and
  conformable-array vectorization in standard Fortran.
- `pbivnorm_recycle` reproduces R's arbitrary-length recycling behavior.
- Infinite bounds and the singular correlations `rho = +/-1` are handled
  explicitly rather than by replacing infinities with `double.xmax`.
- Invalid `abs(rho) > 1` returns IEEE NaN in the scalar/elemental interface;
  `pbivnorm_recycle` additionally reports a nonzero status.
