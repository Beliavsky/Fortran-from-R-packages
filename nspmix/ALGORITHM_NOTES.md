# Algorithm notes

## Constrained Newton mass step

For component density matrix `D`, weights `w`, and current proportions `p`, the
port forms `P=D p`, `S=D/P`, and solves the Wang constrained-Newton quadratic
subproblem under the probability-simplex constraints. Backtracking enforces
monotone log-likelihood increase. This is the computational heart of `hcnm`
and the mass-update step used by `cnm`.

The attached LSEI translation is used through `lsei_solve` with one equality
constraint and lower bounds of zero. This is equivalent to the upstream
partial-NNLS simplex subproblem.

## NPMLE support search

`cnm` evaluates the directional gradient on family-specific support grids,
adds the best positive-gradient support point, refits all masses, removes tiny
components, and repeats until the likelihood/support-gradient stopping rule is
met or `kmax` is reached.

## Semiparametric models

Normal and CVPS use one positive scale parameter. Mlogit uses one fixed-effect
coefficient per covariate. The same analytic component scores as upstream are
implemented. In v0.1.0 the three semiparametric entry points use a common
profile/alternating optimizer instead of reproducing three separate R BFGS
schedules. This is documented rather than hidden because convergence paths can
differ even when the optimized likelihood is the same.

## Numerical safeguards

Mixture likelihoods use log-sum-exp. Poisson and probability-family supports
are clipped only where upstream also protects logs/divisions. All source is
free-form Fortran and avoids relying on short-circuit Boolean evaluation.
