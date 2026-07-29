# VaRES for modern Fortran

This project is a modern Fortran 2018 translation of VaRES 1.0.2, an R
package for value at risk (VaR), expected shortfall (ES), probability density
functions, and cumulative distribution functions for more than one hundred
parametric distributions.

## Scope

The library provides all 448 exported computational procedures from the R
package:

- 110 density functions with names beginning with `d`
- 110 cumulative distribution functions with names beginning with `p`
- 114 quantile/VaR functions with names beginning with `var`
- 114 lower-tail expected-shortfall functions with names beginning with `es`

Four families in the original package provide only VaR and ES. The remaining
110 families provide the complete density/CDF/VaR/ES interface.

All public distribution procedures are `pure elemental`. A scalar call and a
conformable array call therefore use the same API.

## Design

- Modern Fortran 2018 and FPM project layout
- `dp = kind(1.0d0)` throughout
- No external numerical-library dependency
- Self-contained normal, beta, gamma, Student-t, F, lognormal, logistic,
  Cauchy, uniform, and Weibull support functions
- Continued-fraction and series evaluations for incomplete beta and gamma
  functions
- Bracketed inverse-CDF calculations for beta and gamma distributions
- 96-point Gauss-Legendre ES integration with an endpoint transformation
- Original default parameter values and lower-tail conventions
- Original R source and documentation retained under `original/`

## Build and test

With FPM installed:

```text
fpm build
fpm test
fpm run
fpm run --example vector_example
```

The package has no external dependencies.

## Example

```fortran
program example
  use vares, only : dp, pnormal, varnormal, esnormal
  implicit none

  real(dp), parameter :: p = 0.05_dp
  real(dp) :: q

  q = varnormal(p, mu=0.0_dp, sigma=1.0_dp)

  print '(a,f12.6)', 'VaR: ', q
  print '(a,f12.6)', 'ES:  ', esnormal(p)
  print '(a,f12.6)', 'CDF: ', pnormal(q)
end program example
```

For CDF and quantile procedures, the R arguments `log.p` and `lower.tail`
become `log_p` and `lower_tail`. For density procedures, `log` becomes
`log_pdf`.

## Numerical conventions

`var*` returns the lower-tail quantile at probability `p`. `es*` returns

```text
ES(p) = (1/p) integral from 0 to p of VaR(u) du.
```

This is the convention used by the original package. Some heavy-tailed
families do not have a finite ES for every parameter choice; in those cases a
large magnitude, infinity, or NaN can be mathematically appropriate.

## Source corrections

Several original R formulas do not satisfy their own density/CDF/quantile
identities or ignore `log.p` and `lower.tail`. The Fortran port corrects those
cases rather than preserving demonstrably invalid results. See `PORTING.md`
for the complete list and the exact rationale.

## Documentation

- `API.md`: complete generated procedure inventory
- `PORTING.md`: R-to-Fortran mapping and source corrections
- `TESTING.md`: test coverage and compiler commands
- `NOTICE.md`: original authorship, citation, and provenance
- `CHANGELOG.md`: port release notes

## License

The original package declares `GPL (>= 2)`. This derivative work is therefore
licensed under GPL-2.0-or-later. The complete GPL-2.0 and GPL-3.0 texts are in
`LICENSES/`, and the top-level `LICENSE` contains GPL-2.0.
