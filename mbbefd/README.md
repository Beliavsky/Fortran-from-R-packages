# mbbefd-fortran

Modern Fortran 2018/FPM computational port of the R package **mbbefd 0.8.14**.
The original package is by Christophe Dutang and contributors and is distributed
under GPL-2. This port preserves the license and retains the supplied upstream
source tree under `upstream/mbbefd-master/` for provenance.

## Scope

The port covers the computational portions of mbbefd:

- MBBEFD distributions in both `(a,b)` and `(g,b)` parameterizations
- density/mass, CDF, quantile, RNG, exposure curve, moments and total-loss mass
- shifted truncated Pareto and generalized-beta-first-kind distributions
- one-inflated uniform, beta, shifted-Pareto and generalized-beta distributions
- generic one-inflated distribution helpers using Fortran procedure callbacks
- Swiss Re curve parameter helper and `g2a`
- empirical total loss, Theil statistic and empirical exposure curves
- parameter transforms used by fitting
- analytic MBBEFD log-likelihood gradients and `(a,b)` Hessian
- `fit_dr` with MLE and total-loss/moment matching
- `boot_dr` parametric and nonparametric bootstrap

R graphics, S3 printing/summary methods, `eccomp`, and Rcpp glue are omitted.

## Dependencies

The supplied translations are used as real FPM dependencies and are vendored so
this archive is self-contained:

- `fitdistrplus-fortran` 0.1.0
- `actuar` 0.1.0
- `alabama` 0.1.1 (including its supplied `numDeriv` and `roptim` dependencies)

No R runtime is required.

## Build

```sh
fpm build
fpm test
fpm run --example basic
```

The package was validated directly with gfortran 14.2 because FPM was not
installed in the translation environment. The FPM manifest itself was parsed
and the same `src/`, `test/`, `example/`, and dependency trees were compiled.

## Example

```fortran
program demo
  use mbbefd, only : dp, dmbbefd, pmbbefd, qmbbefd
  implicit none

  print *, dmbbefd(0.4_dp, 0.5_dp, 0.3_dp)
  print *, pmbbefd(0.4_dp, 0.5_dp, 0.3_dp)
  print *, qmbbefd(0.25_dp, 0.5_dp, 0.3_dp)
end program demo
```

The R package names `dMBBEFD`, `pMBBEFD`, etc. differ from `dmbbefd`,
`pmbbefd`, etc. only by case. Fortran identifiers are case-insensitive, so the
`(g,b)` versions are named `dmbbefd_gb`, `pmbbefd_gb`, `qmbbefd_gb`, and so on.
See `API_MAP.md` for the complete mapping.

## Numerical notes

MBBEFD distributions contain an atom at 1 in addition to a continuous part.
Consequently the `d*` routines return the probability mass when `x == 1`, in
the same convention as the upstream package.

For arbitrary positive raw moments where no short closed form is used, the port
integrates the survival-function identity with fixed high-order Gauss-Legendre
quadrature. See `PORTING_NOTES.md` for compatibility details and deliberate
fixes to upstream edge cases.
