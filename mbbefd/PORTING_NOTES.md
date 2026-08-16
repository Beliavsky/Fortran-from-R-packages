# Porting notes

## Design

The port translates the numerical/statistical code to Fortran 2018 and keeps
the supplied Fortran translations of `fitdistrplus`, `actuar`, and `alabama`
as explicit local FPM dependencies. Dependency algorithms are not copied into
`mbbefd` source modules.

The public modules are:

- `mbbefd_kinds`
- `mbbefd_math`
- `mbbefd_distributions`
- `mbbefd_fit`
- umbrella module `mbbefd`

## Case-insensitive identifiers

R exports two MBBEFD parameterizations as lower-case `mbbefd` and upper-case
`MBBEFD`. Fortran cannot distinguish identifiers by letter case. The `(a,b)`
API retains the lower-case spelling while the `(g,b)` API uses the suffix
`_gb`.

## Fitting

`fit_dr` supports the upstream `mle` and `tlmme` methods. One-inflated beta,
generalized-beta and shifted-Pareto MLEs use `fitdistrplus-fortran`. The
MBBEFD families have disconnected parameter domains, so both admissible regions
are optimized with `alabama` and the better converged solution is selected.

The supplied `alabama-fortran` callback API does not carry a user context.
Therefore the small MBBEFD/alabama adapter stores the observations in module
state for the duration of a constrained optimization. Calls to this specific
fit path must not be executed concurrently in the same process. Ordinary
`fitdistrplus` fits do not have this restriction.

The supplied `alabama-fortran`/`roptim` tree also emits legacy labeled-DO
compiler warnings, and GNU ld reports an executable-stack request originating
inside `alabama.o`. These diagnostics are in the supplied dependency rather
than the translated mbbefd source. The constrained-fit runtime tests were run
successfully with that dependency as supplied.

## Compatibility fixes

A few upstream implementation details were corrected rather than copied
literally:

1. The generic R `qoifun` scales the input probability before interpreting
   `lower.tail` and `log.p`, and then compares the already-scaled probability
   against `1-p1`. The Fortran implementation first converts the user's input
   to an ordinary lower-tail probability and applies the mathematically correct
   mixture quantile: return 1 for `p >= 1-p1`, otherwise evaluate the base
   quantile at `p/(1-p1)`.
2. `constrmbbefd1jac` in the R source has inconsistent derivative rows. The
   Fortran Jacobian differentiates the actual four constraints
   `(a+1, -a, b-1, a*(1-b))`.
3. `constrmbbefd2` repeats `a` in its second constraint even though its comment
   and Jacobian specify the lower bound on `b`. The Fortran constraint is
   `(a, b, 1-b, a*(1-b))` and uses its exact Jacobian.
4. One-inflated beta `d/p/q/r` accepts noncentrality just as the R API does.
   Exposure curves and moments remain central-beta-only because the upstream
   package likewise does not implement the noncentral versions.
5. The R shifted-truncated-Pareto and MBBEFD moment functions use numerical
   integration for non-unit orders. The Fortran port exposes arbitrary positive
   real raw moments and uses high-order Gauss-Legendre integration.

## EECF

R returns an `eecf` closure with S3 metadata. The Fortran equivalent is a
small derived type containing the sorted sample and its mean. Calling
`obj%evaluate(d)` computes `mean(min(X,d))/mean(X)`.

## Validation

The tests cover fixed density/CDF/quantile reference values, equivalent MBBEFD
parameterizations, generalized beta, shifted Pareto, exposure/moment helpers,
RNG means and atom masses, one-inflated fitting, constrained MBBEFD fitting,
and parametric bootstrap. Both strict runtime-checked and optimized builds are
part of release validation.
