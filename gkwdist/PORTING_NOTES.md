# Porting notes

## Scope

`gkwdist` 1.1.4 exports 50 computational functions plus the imported `%>%`
operator. All 50 computational exports are represented. There is no plotting
algorithm in the core exported API; R documentation examples that create plots,
Rcpp registration code, Armadillo plumbing, and the pipe operator are omitted.

## Unified GKw implementation

All seven distributions are nested cases of the five-parameter generalized
Kumaraswamy distribution. The Fortran implementation evaluates one stable GKw
kernel and obtains the sub-families by fixing parameters:

- BKw: lambda = 1
- KKw: gamma = 1
- EKw: gamma = 1, delta = 0
- McDonald: alpha = beta = 1
- Kw: gamma = 1, delta = 0, lambda = 1
- upstream `beta_`: alpha = beta = lambda = 1, so the ordinary beta shapes are
  `(gamma, delta + 1)`

This removes duplicate source while retaining the same mathematical formulas.

## Analytical derivatives

The upstream C++ source hand-codes a gradient and Hessian for each nested
family. The Fortran translation evaluates the same closed-form negative
log-likelihood with a second-order automatic-differentiation scalar carrying
its value, gradient, and Hessian. This is analytical differentiation, not
finite differencing. It also makes all nested-family derivatives projections of
the same GKw likelihood, preventing formula drift among seven duplicated
implementations.

The supplied `numDeriv-fortran` package is a validation dependency only. The
test suite compares the analytical Fortran gradient and Hessian with Richardson
numerical derivatives for all seven families.

## Numerical functions

R Mathlib/Rcpp dependencies were replaced with standalone Fortran routines for
regularized incomplete beta probabilities, beta quantiles, log-beta, digamma,
trigamma, gamma RNG, and beta RNG. GKw calculations retain the upstream
log-space `log(1-exp(.))` strategy; the AD implementation includes a
differentiated stable `log1mexp` primitive.

Upper-tail log-quantile conversion uses a stable `-expm1(log_p)` equivalent,
avoiding cancellation when an upper-tail probability is tiny.

## Starting values

`gkwgetstartvalues` retains the upstream design: sample moments 1 through 5,
weights `(1, .8, .6, .4, .2)`, numerical theoretical moments, multiple starts,
and Nelder-Mead minimization. The first four deterministic starts and the
family-specific final bounds follow the C++ implementation. Extra random starts
use a deterministic local Park-Miller stream rather than C++ `mt19937`, so
exact start vectors can differ while the algorithm and objective are the same.

## Intentional compatibility corrections

- Invalid likelihood data outside `(0,1)` return positive infinity consistently
  with a negative-log-likelihood minimization objective. One current GKw C++
  branch returns negative infinity despite its own documentation specifying
  positive infinity.
- The translated density does not silently truncate finite observations merely
  because `x**alpha` lies within a square-root-epsilon band of one. This follows
  the package NEWS statement that the near-boundary truncation was removed and
  lets the log-space formula determine the value.

## Array and RNG semantics

Fortran elemental procedures support scalar broadcasting and conformable array
evaluation. Arbitrary R recycling between nonconformable vector lengths is not
performed. Random generation uses Fortran `random_number` plus a
Marsaglia-Tsang gamma generator, so random streams do not match R's RNG bit for
bit.
