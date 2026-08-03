# Porting notes

## R objects

R lists and attributes are represented by the derived types `sample_moments`,
`skew_t_parameters`, and `portfolio_result`. Optional arguments replace R's
named list components and `...` mechanism.

## Co-moment representation

The R package expands compressed PerformanceAnalytics co-moments into matrices
used by SCA algorithms. A full fourth-order object scales as p^4. The Fortran
port instead retains centered observations and computes portfolio moments,
gradients, and Hessians from matrix-vector operations. This is algebraically
equivalent for the sample central moments and makes the default representation
usable for substantially larger portfolios. Full tensors are optional.

## Skew-t dependency

The required generalized hyperbolic skew-t EM/PX-EM estimator was copied and
adapted from the earlier GPL-3.0-only `fitHeavyTail-fortran` translation. Only
the modules required by `fit_mvst` are included.

## Optimization dependencies

The upstream package uses `quadprog`, `ECOSolveR`, `lpSolveAPI`, and `nloptr`.
This project is dependency-free:

- convex simplex QPs use accelerated projected gradient;
- `Q-MVSK`, `MM`, and `DC` retain their SCA/majorization structures;
- `PGD`, `RFPA`, and `SQUAREM` retain projected-gradient and acceleration logic;
- tilting directly maximizes the worst signed moment improvement and projects
  trial portfolios onto the tracking ellipsoid along the reference-portfolio
  ray.

The tilting subproblem solver is therefore a numerical adaptation, not a call
to the exact ECOS/LP/QP subproblems used by R. It preserves the objective,
simplex, and tracking-error constraints.

## Magnitude adjustment

`adjust_magnitude=.true.` rescales each moment so the equal-weight portfolio has
absolute moment one. The sign of odd moments is retained, matching the R code.

## Omitted material

Plotting is not part of the package's exported computational API. R package
metadata, RData datasets, vignettes, and precompiled `.o`/`.so` files are not
part of the Fortran build. Relevant R source, manuals, tests, and metadata are
retained in `original/`.
