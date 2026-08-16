# poweRlaw-fortran

Modern Fortran 2018/FPM translation of the computational core of the R package
`poweRlaw` 1.0.0 (Analysis of Heavy Tailed Distributions).

The port provides the eight fitted distribution families from the R package:

- discrete and continuous power law
- discrete and continuous exponential
- discrete and continuous lognormal
- discrete Poisson
- continuous Weibull

It also translates parameter estimation, KS/reweighted lower-cutoff (`xmin`)
estimation, Vuong distribution comparison, ordinary bootstrap, and the
semi-parametric power-law plausibility bootstrap.

## FPM

```text
fpm build
fpm test
fpm run --example basic_fit
```

The supplied `pracma-fortran` 2.4.6 translation is vendored under
`vendor/pracma-fortran`.  The power-law normalizer uses an Euler-Maclaurin
Hurwitz-zeta implementation for the difficult alpha-near-one region and reuses
`pracma`'s Riemann zeta implementation in its rapidly convergent high-alpha
case.

## Minimal example

```fortran
program demo
   use powerlaw
   implicit none
   type(powerlaw_dist) :: m
   type(estimate_xmin_result) :: fit
   real(dp) :: x(20)
   integer :: i

   do i = 1, size(x)
      x(i) = real(i,dp)
   end do

   m = displ(x)
   fit = estimate_xmin(m)
   print *, fit%xmin, fit%pars, fit%gof
end program
```

## Design

R reference classes are represented by `type(powerlaw_dist)`.  Constructors
such as `displ`, `conpl`, and `conlnorm` select the family.  The object stores
sorted data, `xmin`, and the fitted parameters and exposes type-bound PDF/CDF,
log-likelihood and random-generation operations.  Fitting and bootstrap calls
return typed result structures rather than R lists/data frames.

Bootstrapping is serial in v0.1.0; the R package's `parallel` orchestration is
not part of the numerical translation.  Plotting, S4/reference-class display
methods and other R presentation machinery are intentionally omitted.

See `API_MAP.md` and `PORTING_NOTES.md` for exact coverage and numerical notes.
