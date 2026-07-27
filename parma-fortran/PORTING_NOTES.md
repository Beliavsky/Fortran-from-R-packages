# Porting notes

## Language mapping

R vectors and matrices map to `real(dp)` arrays. Heterogeneous S4 objects map to
`parma_spec`, `parma_port`, and solver result derived types. R functions passed
as optimization objectives map to the explicit `objective_callback` abstract
interface. All public Fortran procedures have explicit interfaces and all
sources use `implicit none`.

## Numerical choices

The original package delegates most optimization to `nloptr`, `quadprog`,
`Rglpk`, and a bundled C SOCP implementation. To keep the FPM project
self-contained, the translation supplies native algorithms:

- Full-covariance CMA-ES for general nonlinear and nonsmooth portfolios.
- Projected-gradient convex QP.
- Canonical simplex LP.
- Exact small binary enumeration.
- Newton logarithmic-barrier SOCP.

Results need not be bit-for-bit identical to a particular external R solver,
especially for nonunique optima. The objective definitions and constraints are
preserved at the mathematical level.

## Source inconsistencies handled explicitly

1. The upstream `sentropy` and `centropy` checks use expressions such as
   `if (any(w) < 0)`, which test a logical result rather than the weights. The
   Fortran routines validate each weight directly.
2. Public `fun.lpm` computes `max(return-threshold,0)`, while the NLP optimizer
   computes `max(threshold-return,0)`. The Fortran default is the conventional
   downside partial moment. `legacy=.true.` and `spec%lpm_legacy=.true.` retain
   the public R expression when compatibility is required.
3. The documented threshold value `999` now consistently centers portfolio
   returns before LPM/UPM evaluation.
4. CDaR is implemented from the standard Rockafellar-Uryasev tail-average
   expression over drawdowns. This avoids fragile index/reversal operations in
   the R helper.
5. The R matrix lag helper refers to an undefined variable `i` in its scalar-lag
   branch. The Fortran routine uses the requested lag directly.
6. Benchmark-relative covariance risk is evaluated directly as
   `sigma_b^2 + w'Sw - 2 w'cov(asset,b)`.

## Array ordering

Fortran is column-major, matching R's underlying matrix order. `reshape` calls
in tests and examples therefore list each column contiguously.

## Randomness

`seed_rng` expands one integer into the processor's full `random_seed` vector.
CMA-ES and `simweights` are reproducible for a fixed compiler/runtime and seed,
but streams are not intended to duplicate R's Mersenne Twister.

## SOCP

`socp_solve` accepts the same mathematical form documented by the original
`Socp` function. It returns the primal solution and residual information, not
the original C solver's dual vectors or iteration-history workspace.

## License and provenance

The supplied source archive checksum is in
`original/SOURCE_ARCHIVE_SHA256.txt`. The complete original tree is retained
under `original/parma-1.7`. No original copyright notice was removed.
