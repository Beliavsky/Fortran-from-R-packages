# Porting notes

## Source package

- Package: VaRES
- Version: 1.0.2
- Publication date: 2023-04-22
- Original language: R
- Original license: GPL (>= 2)
- Fortran license expression: GPL-2.0-or-later

The original source and manual pages are retained under `original/` for
provenance and formula comparison.

## API mapping

The original exported names are retained. Fortran is case-insensitive, and the
source uses lower-case spelling consistently.

| R argument | Fortran argument |
| --- | --- |
| `log` | `log_pdf` |
| `log.p` | `log_p` |
| `lower.tail` | `lower_tail` |

All scalar R formulas become `pure elemental` Fortran functions. This supports
scalars and conformable arrays without separate vector wrappers. R vector
recycling is not reproduced; arrays should be conformable.

The original ES procedures repeatedly call R's `integrate()` over the
corresponding quantile function. The Fortran implementation applies the same
definition using fixed 96-point Gauss-Legendre quadrature after the transform
`u = p*s**3`. The transform improves accuracy near an infinite or singular
lower endpoint while avoiding probabilities so small that `1-u` rounds to one
in double precision.

## Corrected source inconsistencies

The following changes restore mathematical identities or documented flag
behavior. They are intentionally different from the executable R source.

1. `pinvexpexp`: removed an extra factor `a` from the lower-tail log CDF.
2. `pnakagami`: now honors all `log_p` and `lower_tail` combinations.
3. `varclg`: corrected the lower-branch location sign so `pclg(varclg(p))=p`.
4. `plaplace`: corrected four log-CDF conditions that repeated
   `log.p==FALSE` and overwrote the ordinary CDF.
5. `varloglaplace`: corrected the upper-branch normalization from `a` to `b`.
6. `pasypower`: replaced gamma densities with gamma CDFs and corrected the
   right-side signs and complements.
7. `pgenlogis3`: now passes `log_p` and `lower_tail` to the beta CDF.
8. `pgenlogis4`: now passes `log_p` and `lower_tail` to the beta CDF.
9. `pstacygamma`: now passes `lower_tail` to the gamma CDF.
10. `pkumloglogis` and `varkumloglogis`: corrected the CDF orientation and
    quantile algebra so the density is nonnegative and CDF/quantile invert.
11. `pburr7`: removed the erroneous minus sign from the lower-tail log CDF.
12. `vardagum`: changed `1-p**(-1/c)` to `p**(-1/c)-1`, avoiding a negative
    fractional-power base and restoring inversion.
13. `dkumweibull`, `pkumweibull`, and `varkumweibull`: made the implementation
    consistent with the Kumaraswamy-Weibull density and formulas documented in
    the package manual.
14. `dexpweibull` and `pexpweibull`: use the Weibull CDF rather than its
    survival function, matching the documented exponentiated-Weibull density,
    CDF, and quantile.

The targeted correction tests check CDF/quantile inversion, density/CDF
numerical derivatives, tail complements, and log-probability behavior.

## Special functions

The port is self-contained. The incomplete beta and incomplete gamma functions
use standard continued-fraction and series algorithms. Quantiles use monotone
bracketing and bisection where no direct approximation is provided. This makes
the package portable and avoids requiring a particular BLAS, LAPACK, or
special-function library.

The exact last bits can differ from R because R delegates several functions to
its platform math library. Tests use tolerances appropriate to double
precision and the implemented algorithms.

## Error behavior

The original formulas perform little explicit parameter validation. The port
largely retains that behavior: mathematically invalid parameter combinations
can produce IEEE NaN or infinity. Invalid ES probabilities outside `(0,1]`
return quiet NaN.

## Omitted features

No computational export was omitted. R-only package machinery, help-system
execution, vector recycling, and R-specific NA attributes are not part of the
Fortran API. The original manual pages remain available under `original/man/`.
