# Porting notes

## Objectives and callbacks

An R call such as:

```r
maxLik(logLik, grad=gradient, hess=hessian, start=start, method="BFGS")
```

becomes:

```fortran
call initialize_problem(problem, size(start), objective)
problem%gradient => gradient
problem%hessian => hessian
call max_lik(problem, start, result, 'bfgs')
```

The Fortran objective always returns one scalar total objective. If an R
objective returned one contribution per observation, sum those contributions in
the objective callback and return their individual score rows through the
`scores` callback.

## Parameter transformations

Named R vectors are positional Fortran arrays. Keep parameter positions in a
small module or define integer named constants. Positive parameters should
usually be optimized on a log scale, as in the normal-MLE example.

## Fixed parameters

Use `set_fixed(problem, [indices])`. Solvers retain the starting values at fixed
positions. Final covariance rows and columns for fixed positions are zero.

## Constraints

The upstream convention is retained:

```text
A * theta + b = 0
A * theta + b >= 0
```

The Fortran implementation uses successively larger quadratic penalties.
Always inspect `result%constraint_violation`. Increase
`constraint_max_outer`, `constraint_rho_factor`, or relax `constraint_tol` if a
problem is badly scaled.

## Derivatives

Analytic derivatives are recommended for speed and stability. When omitted,
central finite differences are used by default. Run `compare_derivatives`
before fitting a new model. BHHH and stochastic methods require observation
scores even if a full analytic gradient is also supplied.

## Scaling

Newton and BHHH solve linear systems without external BLAS/LAPACK. Rescale
parameters and objective components when magnitudes differ greatly. Use
`condition_number` or `progressive_condition_numbers` to diagnose unstable
parameterizations.

## Random methods

Simulated annealing, SGA, and Adam use a self-contained seeded generator.
Set `control%random_seed` for reproducible runs. Stochastic methods return the
best full-objective parameter vector observed, not necessarily the final
mini-batch iterate.
