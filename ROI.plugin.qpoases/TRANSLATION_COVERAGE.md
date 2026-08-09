# Translation coverage

## Translated

### ROI.plugin.qpoases layer

* quadratic objective convention `0.5*x^T*Q*x + L^T*x`
* minimization/maximization sign conversion
* ROI default nonnegative variable bounds
* `>=`, `<=`, and `==` linear constraint conversion
* automatic ZERO / IDENTITY / UNKNOWN Hessian classification
* control/options structure
* primal solution, dual solution, objective and status output

### qpOASES-facing solver surface

* QProblem-like general QP solves
* QProblemB-like simply bounded QPs
* SQProblem-like persistent matrix/gradient hotstarts
* Hessian type constants 0--6
* lower/upper variable bounds
* lower/upper/equality general constraints
* feasible active-set iterations
* working-set add/drop logic
* KKT multipliers
* Phase-I feasibility QP
* warm reuse of the previous active set
* objective/primal/dual getters
* active/free/fixed/equality count diagnostics
* infeasible status
* zero-Hessian LP recession-direction unbounded detection
* dense linear algebra without external BLAS/LAPACK

## Deliberate architectural differences

The original bundled qpOASES 3.2 implementation is a highly optimized online
active-set code with specialized TQ/QR/Cholesky updates and a detailed
homotopy path.  v0.1.0 preserves the active-set solver model and public
computational behavior but recomputes dense KKT systems instead of
incrementally updating qpOASES' matrix factorizations.

Consequently the following qpOASES internals are not yet one-for-one ports:

* ramping and far-bound homotopy mechanics
* flipping-bound machinery
* drift correction
* full LI/NZC tests and constraint-dropping priorities
* inertia correction
* incremental TQ/QR/Cholesky factorization updates
* sparse matrix / Schur-complement backends
* exact qpOASES cycling-resolution heuristics
* guessed `Bounds` / `Constraints` objects
* OQP benchmark file reader/runner
* `SolutionAnalysis`
* exact CPU-time interruption behavior

All fields of the R package's default option list are retained in
`qpoases_options`; fields belonging only to the omitted advanced internals are
currently informational.

The solver is intended for convex QPs.  The obvious zero-Hessian LP
unboundedness case is detected through a recession-direction working-set
solve.  General indefinite-QP global optimization is not implemented; this
matches the fact that qpOASES is primarily a convex QP solver despite exposing
Hessian classification constants.

## R-only code omitted

* ROI solver registration/signatures
* R S4/ROI model objects
* Rcpp external pointers and registration
* status database registration
* console printing
* `dry_run` call objects

The complete original source tree is preserved under
`original/ROI.plugin.qpoases-master/`.
