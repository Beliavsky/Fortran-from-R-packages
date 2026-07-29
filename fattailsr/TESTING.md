# Testing

Validation was performed with GNU Fortran 14.2.0 using:

```text
-std=f2008 -Wall -Wextra -Werror -Wimplicit-interface
-fcheck=all -ffpe-trap=invalid,zero,overflow
```

The automated test program covers:

- logit/inverse-logit, incomplete beta, and normal quantile identities;
- all principal parameter conversions;
- K4 quantiles cross-checked against independently calculated constants;
- CDF/quantile and value/logit inversion;
- density versus a finite-difference CDF derivative;
- reciprocal density/quantile derivatives;
- left/right tail means and signed/unsigned expected shortfall;
- theoretical raw and central moments through order four;
- exact synthetic recovery by the five-, seven-, and eleven-quantile
  estimators;
- exact synthetic recovery by the K4 fitter;
- price return, missing-value, and elevation transformations;
- finite, nonconstant random generation.

Result:

```text
All FatTailsR tests passed.
```

The environment used to prepare this archive did not contain the `fpm`
executable. The standard FPM directory structure and manifest were supplied,
and all FPM targets were compiled and run directly with `gfortran`.
