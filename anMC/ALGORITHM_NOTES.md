# Algorithm notes

## Orthant decomposition

`proba_max` and `proba_min` follow the upstream decomposition

`p = p_q + R_q (1 - p_q)`

where `p_q` is computed on a selected active subset and `R_q` is a conditional
remainder probability estimated by ordinary MC or ANMC.  Conditional Gaussian
means/covariances are formed with the same Schur-complement formulas as the R
code.

The upstream variance combination is preserved:

`(1-R_q)^2 Var(p_q) + Var(R_q)(1-p_q)^2 + Var(R_q)Var(p_q)`

with `Var(p_q) = (pmvnorm_error/3)^2` in the bias-corrected branch.

## ANMC

The timing-based calibration of outer cost, conditional-simulation cost, and
indicator-evaluation cost is retained.  The initial `(n0,m0)` experiment is
used to estimate the conditional variance ratio and derive `mStar`; the final
outer count is then obtained from the remaining computational budget.

The published nested-MC variance estimator used by the upstream code is
preserved.

## Active dimensions

All six upstream heuristics are implemented:

0. equally spaced indices;
1. weighted sampling by `pn`;
2. weighted sampling by `pn(1-pn)`;
3. `pn` weighted by accumulated design-space distance;
4. `pn(1-pn)` weighted by accumulated design-space distance;
5. uniform sampling without replacement.

`select_q_dims` increases `q` by `min(10, ceil(0.01*n))` until successive active
probability estimates differ by no more than the integration error or the
specified upper limit is reached.

## Native Gaussian probability kernel

`anmc_math::pmvnorm` standardizes the rectangle to a correlation problem,
orders dimensions by marginal interval probability, factors the correlation
matrix, and evaluates the sequential Genz conditional transform with randomized
Halton points.  Random shifts use a private integer generator, so probability
integration does not reset or consume the Monte Carlo simulation RNG stream.

## Intentional corrections to upstream code

1. In `ANMC_Gauss`, upstream code sets `indM=min(m0,mStar)` and then tests
   `if(indM>m0)` before generating missing conditional draws.  That condition is
   impossible.  The port uses the intended `if(mStar>m0)` condition.
2. `conservativeEstimate` allocates `productPn` with fixed length 10000 and can
   index beyond both that vector and the sorted probabilities if every marginal
   probability exceeds `alpha`.  The port uses bounded/dynamic logic.
3. `conservativeEstimate` documents algorithm labels `GANMC` and `GMC`, while
   upstream `ProbaMax`/`ProbaMin` only test for the exact string `ANMC`; therefore
   the default `GANMC` falls through to ordinary MC.  The port maps both `ANMC`
   and `GANMC` to ANMC, and `MC`/`GMC` to ordinary MC.
4. One large-set branch passes `E[1:NextEval]` rather than `E[1:NextEval,]`,
   flattening a multidimensional R design.  The Fortran port always preserves
   the design matrix shape.
5. Active-covariance reselection is bounded to 100 attempts instead of allowing
   an unbounded loop on a persistently singular problem.
6. Distance-weighted active selection falls back safely when all accumulated
   distance weights collapse to zero rather than propagating NaNs into weighted
   sampling.
7. Rejection-sampling batch sizes are explicitly capped before integer
   conversion, preventing overflow or accidental enormous allocations at very
   small estimated acceptance probabilities.
8. Indicator timing in ANMC measures the actual vector `maxval`/`minval` work;
   the R code times `gg(1)`, a scalar operation, even when the conditional vector
   is high dimensional.

These corrections do not change the defining probability decomposition or MC
estimators.
