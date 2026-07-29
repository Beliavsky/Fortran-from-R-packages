# Testing

The translation was tested with GNU Fortran 14.2.0 using:

```text
-std=f2018
-Wall -Wextra -Werror
-fcheck=all
-fbacktrace
-ffpe-trap=invalid,zero,overflow
```

`test/test_dowd.f90` covers:

- normal and Student-t CDF/quantile behavior
- independently calculated normal, Student-t, lognormal VaR and ES values
- historical and Cornish-Fisher ordering
- Gumbel, Frechet, GPD, Hill, and Pickands calculations
- KDE, bootstrap, and Box-Cox routines
- portfolio VaR/ES, hotspots, eigendecomposition, and full-factor PCA identity
- coverage and independence backtests and goodness-of-fit statistics
- Black-Scholes call/put reference prices
- American versus European put price ordering and American VaR/ES ordering
- product, Gaussian, and Gumbel copula identities

The strict direct build and test command completed with:

```text
All Dowd tests passed.
```

The demonstration and option-risk example also compiled and ran successfully.
An `fpm` executable was not installed in the build environment, so the standard
FPM project structure and manifest were compiled directly with `gfortran`.
