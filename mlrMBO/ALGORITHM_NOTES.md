# Algorithm notes

## Surrogate model

Upstream mlrMBO can use many regression learners through `mlr`. The Fortran
port uses the translated DiceKriging Gaussian-process model as its native
surrogate. One independent model is fitted per objective. Universal-kriging
prediction supplies both posterior mean and standard error. The covariance
family, multi-start controls and covariance re-estimation behavior are
available through `mbo_control`.

A tiny fixed nugget (`1e-10`) is used during fitting to make repeated or
nearly repeated design points numerically safer. Proposed-point filtering is
also enabled by default.

## Infill criteria

The formulas follow `R/infill_crits.R`. As in mlrMBO, all criteria are
converted to an internal minimization problem. Consequently EI and standard
error are returned with a minus sign. CB is the lower confidence bound for a
minimized objective (or negative upper confidence bound for maximization).

AEI and EQI require a pure-noise estimate. When no explicit model nugget is
requested, the Fortran port estimates this from training residuals. This is
the closest standalone analogue of mlrMBO's fallback to residual-variance
estimation through its learner framework.

## Focus search

The native focus optimizer follows `infillOptFocus.R`: independent restarts,
Latin-hypercube candidate sets, and repeated halving of the numerical search
range around the best local candidate. Continuous and integer ranges shrink;
categorical level sets are not pruned during focus iterations.

## Multi-objective calculations

The SMS and additive epsilon formulas are direct translations of the small C
kernels in `src/infill.c`. Dominated hypervolume is independently implemented
with recursive slicing of the union of axis-aligned boxes. It accepts any
objective dimension and filters dominated/out-of-reference points first.

DIB converts every objective to minimization, forms confidence-bound vectors,
then minimizes either the SMS or epsilon indicator. For batch DIB the
predicted CB vector of each proposal is inserted into a temporary front before
selecting the next point, following the intended upstream sequential logic.

ParEGO uses range normalization and the augmented Tchebycheff scalarization
`max(lambda*y) + rho*sum(lambda*y)`. Weight vectors are sampled from the same
stars-and-bars grid parameterized by `parego_s`, rather than from a continuous
Dirichlet distribution.

Upstream MSPOT calls an external NSGA-II implementation. The Fortran port
replaces that search engine with a large LHS candidate pool, nondominated
surrogate confidence-bound filtering, and removal-hypervolume selection. The
surrogate objectives are the same class of quantities, but exact candidate
trajectories are therefore not expected to match R.

## Batch proposals

Constant liar repeatedly optimizes the selected infill criterion and updates
a temporary GP with the chosen lie. Parallel-CB draws independent Exp(1)
lambdas and optimizes each confidence bound. Execution is sequential in this
standalone library; the resulting proposal definition is independent of
whether evaluations are later run concurrently.

MOI-MBO is deliberately deferred because the upstream method is an explicit
`emoa` (mu+1) evolutionary algorithm using SBX and polynomial mutation. The
reserved selector raises an error rather than pretending another diversity
heuristic is the same algorithm.
