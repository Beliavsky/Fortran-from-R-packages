# Porting notes

## Source and license

The source package is `sde` 2.0.21 by Stefano Maria Iacus, published under
`GPL (>= 2)`. This derivative work therefore uses `GPL-2.0-or-later`. The
original `DESCRIPTION`, `NAMESPACE`, and `NEWS` files are retained in
`upstream/`. Each translated source file has an SPDX license identifier.

## From R expressions to Fortran procedures

The R implementation accepts expressions and evaluates them in environments.
The Fortran implementation replaces those expressions with explicit procedure
interfaces. This gives compile-time type checking and removes dependence on an
interpreter.

The main coefficient interface is:

```fortran
pure function sde_coefficient(t, x, theta) result(value)
   real(dp), intent(in) :: t, x, theta(:)
   real(dp) :: value
end function sde_coefficient
```

Time-independent information and exact-acceptance routines use
`state_function(x, theta)`. Other interfaces are declared in
`src/sde_interfaces.f90`.

Symbolic differentiation is not attempted. Derivative callbacks such as
`drift_x`, `drift_xx`, `diffusion_x`, and `diffusion_xx` must be supplied where
required.

## Data and time conventions

R `ts` metadata is replaced by explicit `dt`, `t0`, and `t_end` arguments.
Paths are stored with time in the first dimension:

```text
path(step_index, path_index)
```

The initial condition is row 1, so a simulation with `n_steps` returns
`n_steps + 1` rows.

## Parameter conventions

OU:

```text
dX = (theta(1) - theta(2)*X) dt + theta(3) dW
```

GBM/Black-Scholes:

```text
dX = theta(1)*X dt + theta(2)*X dW
```

CIR:

```text
dX = (theta(1) - theta(2)*X) dt + theta(3)*sqrt(X) dW
```

## Optimization

R's `optim` and L-BFGS-B are replaced by an internal bounded Nelder-Mead
routine. Result types retain the estimate, objective, iteration/evaluation
counts, status, and convergence flag. GMM also returns its final weighting
matrix, moment mean, and numerical Hessian.

This optimizer is convenient and dependency-free, but it is not a bitwise or
algorithmic reproduction of `optim`. Difficult likelihoods may benefit from a
specialized external optimizer in downstream applications.

## Numerical support

The standalone library includes:

- Acklam inverse-normal approximation with a Halley refinement.
- Gamma CDF and quantile routines based on series/continued fractions and
  bracketing.
- Noncentral chi-square density/CDF/quantile routines based on Poisson mixtures.
- Marsaglia-Tsang gamma generation and transformed-rejection Poisson generation.
- Adaptive Simpson quadrature.
- Pivoted Gaussian elimination, matrix inversion, and finite-difference
  Hessians.

## Documented source corrections and clarified behavior

The translation preserves the intended algorithm while correcting or making
explicit several inconsistencies found in the supplied upstream source:

1. The generic simulator's second drift derivative referenced `drift.x` rather
   than `drift.xx` in one R closure. The Fortran second-order methods use the
   explicit `drift_xx` callback.
2. The R wrapper and C argument naming for Euler predictor-corrector `alpha`
   and `eta` are inconsistent. The Fortran routine gives each control an
   explicit named argument and applies it according to the displayed formula.
3. The upstream stationary CIR random generator uses gamma parameters that do
   not match its stationary density/CDF/quantile routines. The Fortran random
   generator uses the distribution consistent with those routines:
   shape `2*theta(1)/theta(3)^2` and scale
   `theta(3)^2/(2*theta(2))`.
4. The Pedersen simulation loops in the upstream C code have ambiguous
   off-by-one indexing. The Fortran implementation performs
   `n_substeps - 1` simulated interior Euler steps and evaluates the final
   Euler Gaussian density over the remaining substep.
5. The exact acceptance implementation requires the user to supply the
   transformed-process endpoint sampler and the complete `psi` function. This
   avoids relying on the incomplete default symbolic expression in the R
   wrapper.
6. The upstream Hermite C formula accepts a sixth derivative/moment input but
   does not use it in the final coefficient. The translated formula preserves
   that behavior for compatibility.
7. The upstream `linear.mart.ef` collapses contributions into one scalar even
   for multi-parameter models. The Fortran API retains one equation per
   parameter, which is the natural estimating-equation representation.
8. Coefficients in stepwise schemes are evaluated at the current grid time,
   with next-time evaluation only where the predictor-corrector formula calls
   for it explicitly.
9. `MOdist` returns an R `dist` object containing one triangle. The Fortran
   routine returns the corresponding full symmetric distance matrix and can
   optionally return all estimated operators.

These choices are covered by tests or stated explicitly in the API.

## Reproducibility

Call `seed_rng` before stochastic routines. It seeds the compiler intrinsic
random generator and resets the cached Box-Muller normal draw. Exact sequences
may differ across compilers because the Fortran standard does not prescribe a
specific intrinsic generator.

## Validation status

The included test programs compile with strict Fortran 2018 checks and cover
all major module families. They verify analytical distribution round trips,
path invariants, likelihood finiteness, estimating-equation roots, GMM output,
nonparametric output, changepoint recovery, and B-spline/Markov-distance
properties. Monte Carlo methods are tested with fixed seeds and tolerant
statistical checks rather than exact sample paths.
