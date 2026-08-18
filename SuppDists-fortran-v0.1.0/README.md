# SuppDists-fortran

Modern Fortran translation of the computational core of the R package
**SuppDists 1.1-9.9** (Bob Wheeler; Thorsten Pohlert), packaged for FPM.

## Included distribution families

- inverse Gaussian (Wald): density, CDF, quantile, RNG, moments/statistics
- Kruskal-Wallis statistic and normal-score analogue
- Kendall's tau
- Friedman chi-squared and Spearman's rho
- maximum F ratio
- Pearson sample correlation coefficient
- Johnson SN, SL, SU and SB systems, including quantile and moment fitting
- generalized hypergeometric distributions, including Kemp-Kemp type detection
- expected normal order statistics (AS 177)
- sample mean/SD/skewness/excess-kurtosis helper

The umbrella module is `suppdists`. Public procedure names use lower-case
Fortran spellings such as `dinvgauss`, `pkruskalwallis`, `qjohnson`, and
`dghyper`. Fortran is case-insensitive.

## Build

With FPM:

```text
fpm test
fpm run --example demo_suppdists
```

The translation is standard Fortran 2018 and does not require external
libraries.

## Example

```fortran
use suppdists
real(dp) :: p
p = pinvgauss(1.5_dp, 2.0_dp, 3.0_dp)
print *, p
```

See `TRANSLATION_NOTES.md` for provenance, numerical-method details, and
intentional differences from the R interface.
