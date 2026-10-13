# Numerical API and coverage details

All procedures are public through `network_risk_measures`. `dp` is imported from
`r_kinds` in sibling rfortran-core. Arrays are dense, column-major Fortran arrays.
Callers allocate ordinary output arrays to the documented shape. Do not alias
an input with an output. Optional arguments are best passed by keyword.

## Status codes

| Constant | Value | Meaning |
|---|---:|---|
| nr_ok | 0 | Valid and converged result |
| nr_shape | 1 | Invalid array dimensions or empty network |
| nr_domain | 2 | Invalid value, option, tolerance or control |
| nr_iterations | 3 | Iteration cap reached; last valid approximation returned |
| nr_infeasible | 4 | Incompatible totals, zero-diagonal infeasibility, or no feasible stochastic move |
| nr_random_exhausted | 5 | Supplied uniform stream exhausted; last accepted matrix returned |
| nr_numeric | 6 | Floating-point range failure |

Invalid ordinary-array outputs are initialized to NaN. Reconstruction's last
iterate is retained for iteration/RNG/numerical failures after initialization.
Invalid contagion inputs leave allocatable result components unallocated.
An overflow during contagion leaves partial scenario results and diagnostics.
The unchanged sibling exponential can itself stop on a failed internal LAPACK
solve; such dependency failures cannot be intercepted by this wrapper.

## Risk matrices and criticality

`call risk_matrix(exposures, result, status, buffer, binary, exposure_type, returns)`

- exposures/result: n-by-n; n>=1, finite nonnegative inputs.
- buffer: optional length-n capital vector; required for assets/liabilities;
  each element must be at least sqrt(machine epsilon). Provided buffers are
  validated even for already-normalized risk matrices, as upstream does.
- binary: optional false by default, applies the threshold value >=1.
- exposure_type: assets (default), liabilities, impact or vulnerability.
- returns: impact (default) or vulnerability.

`impact_matrix` and `vulnerability_matrix` have the same arguments except
returns, fixed according to their name. `criticality` also omits returns, has
length-n output, and sums each row of the impact matrix after capping at one.
All are pure.

## Communicability and measures

`call communicability_matrix(x, result, status, terms)` computes exp(x).
Both matrices are n-by-n, finite signed x is accepted. terms=-1 (default)
selects the unchanged sibling exponential. A nonnegative integer terms gives
I + x + ... + x**terms/terms!. terms=0 returns I rather than the source R loop
artifact; high degrees use the equivalent scaled-term recurrence.

`call impact_susceptibility(exposures, result, status, buffer, weights, terms, binary, exposure_type)`

Arguments after status are optional; risk-matrix conventions apply except
binary defaults to true. result is length n. Remove the communicability diagonal
and sum each column. Without weights divide by that impact column's incoming
strength (zero incoming strength gives zero). With length-n finite nonnegative
weights, multiply those column sums elementwise by weights without division.
Weights are not normalized for these measures.

`impact_fluidity` uses the same arguments and scalar output equal to the mean
susceptibility. `impact_diffusion` uses the same arguments and n-by-3 output:
start, intermediate and total. For each vertex remove its outgoing impact row,
recompute communicability, zero the diagonal and clamp numerical negative
differences to zero. Without weights divide changed row sums by the vertex's
original outgoing strength; with weights compute difference times the weight
vector. A vertex with no outgoing strength contributes three zeros.

## Contagion

`call contagion(exposures, result, status, buffer, shocks, weights, method, exposure_type, single_hit, alpha, keep_history, abs_tol, max_it)`

- result: `type(contagion_result)`, allocated internally.
- buffer: optional length n, default ones, used in exposure normalization and
  monetary loss summaries, even with impact inputs.
- shocks: optional n-by-number_of_scenarios finite nonnegative stress, default
  identity (one full initial default per vertex). Stress above one is accepted.
- weights: optional length-n finite nonnegative weights with positive total;
  default equal weights, normalized stably to sum one.
- method: debtrank (default) or threshold; no arbitrary R environment callbacks.
- single_hit: optional false, original DebtRank one-transmission behavior;
  ignored by threshold.
- alpha: optional threshold trigger (0,1], default one.
- keep_history: optional false, allocate history when true.
- abs_tol: optional positive epsilon**0.2, Euclidean stress-increment norm.
- max_it: optional min(1000,10*n), must be positive.

The result contains initial/final(n,scenarios), summary(5,scenarios),
iterations(scenarios), converged(scenarios), and optional
history(n,max_it,scenarios). Read history only through each scenario's iteration
count; later entries are NaN. Summary rows are original weighted capped stress,
additional weighted capped stress, original monetary losses, additional monetary
losses, and additional defaults. Financial units follow buffer for losses.

`call summarize_contagion(initial, final, weights, buffer, summary, status, threshold, cap)`
is the pure numeric equivalent of the registered, unexported summary.contagion
method. threshold defaults to one; cap defaults to one. Additional defaults
count initially subthreshold nodes with final stress > threshold-1e-12.
Losses use raw stress and are deliberately uncapped. Summary weights are
normalized here; stress weights in susceptibility/fluidity/diffusion are not.

The propagation engine accumulates transmitted fractions, computes new stress
as previous + flow*impact, and caps only the next flow at one. Stop is strict
increment_norm < abs_tol. At least one step is taken even for zero shocks,
matching the source repeat loop. Omitted buffer/weights are initialized safely;
the source tries to evaluate nrow(v) before v exists and may error or use an
unrelated parent-environment v. This bug is not preserved.
Nonconverged scenarios retain final/history and
set nr_iterations; the caller must inspect this status.

## Reconstruction

`call max_ent(rowsums, colsums, result, status, max_it, abs_tol, iterations, residual)`

Marginals are finite nonnegative length n with matching totals; result is n-by-n.
max_it defaults 100000, abs_tol defaults 0.001. iterations/residual are optional
output diagnostics. Residual is sum(abs(row discrepancy)) plus sum(abs(column
discrepancy)). The off-diagonal ones prior is alternately rescaled by rows and
columns; stop is residual < abs_tol. Zero-total inputs return zeros immediately.
No R names, verbose output, recycling or rectangular matrices are supported.

`call min_dens(rowsums, colsums, result, status, cost, lambda, k, alpha, delta, theta, remove_prob, max_it, abs_tol, seed, uniforms, iterations, residual, draws_used)`

Defaults: cost=1, lambda=1, k=100, alpha=delta=1/sum(rowsums), theta=1,
remove_prob=0.01, max_it=100000, abs_tol=0.001, seed=1. Cost, alpha, delta,
theta must be finite nonnegative; lambda/remove_prob in [0,1]; k>=0. Scalars
alpha/delta deliberately avoid R's vector recycling. The original objective is
-cost*positive_link_count - alpha**2*sum(row_deficits) - delta**2*sum(col_deficits).
It is linear in deficits, despite the squared coefficients.

Candidate non-diagonal positions follow Fortran/R column-major order. The prior
weight is max(row_deficit/column_deficit, inverse), with zero-deficit pairs given
zero probability. Log scaling avoids overflow in those ratios. In the first k
iterations assign lambda*min(deficits), then use full allocation. Metropolis
acceptance and random link removal follow the upstream update. The stop test is
sum(abs(row_deficits)) < abs_tol, exactly as upstream; the reported residual
also includes column discrepancies. This heuristic can fail to converge and
makes no global sparsity-optimum claim.

Each loop consumes a removal-decision uniform, then a categorical uniform, then
an acceptance uniform only for add proposals, including favorable proposals.
Optional uniforms(:) must contain finite draws in [0,1), consumed without reuse.
Exhaustion returns nr_random_exhausted. Otherwise a private Park-Miller generator
uses modulus 2147483647 and multiplier 16807 with 64-bit integer products; seed
is mapped into 1:2147483646. Identical controls/seeds give reproducible results
without altering intrinsic/global RNG state. This does not reproduce R's
Mersenne-Twister stream, probability sorting, or sample.int algorithm.

Upstream sample(v,1) incorrectly samples 1:v when occupied contains one integer;
this translation always chooses an actual occupied position. Zero totals,
negative values, unbalanced marginals and impossible zero diagonals are checked
before iteration; upstream may error, silently return bad margins, or get stuck.
Feasibility uses each row_i+column_i <= total with numerical slack max(abs_tol/4,
32*epsilon*total), also used for equal totals. Boundary feasible matrices may
converge slowly. This validation is not a proof of convergence for the heuristic.

`matrix_estimation` dispatches to these exact public routines. Arguments start
rowsums, colsums, result, status, method, max_it, abs_tol, followed by the
minimum-density options and diagnostics. method is me (default) or md.
All reconstruction routines are pure, including the local stateful RNG helper.
