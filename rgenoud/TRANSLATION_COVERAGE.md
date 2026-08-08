# Translation coverage

## Translated

The active numerical path behind current `genoud()` is translated to standalone
modern Fortran:

- population generation and rank-based reproduction;
- operator-weight normalization and population-size adjustment;
- all evolutionary operators used by floating-point optimization;
- integer variants of the evolutionary operators;
- boundary-domain handling;
- starting values and elitist retention;
- scalar minimization/maximization;
- lexical/vector-valued optimization with customizable ordering;
- MemoryMatrix-style duplicate caching;
- BFGS refinement of the best individual;
- operator-9 local-minimum crossover and `P9mix`;
- numerical/analytical gradient handling;
- gradient-based stopping logic;
- numerical Hessian output;
- sample/population moments.

## Replaced by native Fortran equivalents

Several upstream helper files implement operations that Fortran already
provides directly. Matrix multiplication, matrix-vector multiplication,
transpose/copy, and manual C allocation helpers are represented by Fortran
arrays, `matmul`, `transpose`, and allocatable storage instead of one-for-one
wrapper routines.

`change_order.cpp` and its associated old matrix utilities are retained in the
original source tree but are not called by the current upstream `genoud()`
execution path and are therefore not separately exposed.

## Deliberate differences

1. **R callback/glue layer**
   `.Call`, `SEXP`, environments, registration, R printing, and R exception
   handling are omitted. Objective and derivative routines are normal Fortran
   procedures with explicit interfaces.

2. **`stats::optim` dependency**
   Upstream BFGS/P9 delegates to R's `optim`. The Fortran package supplies a
   self-contained inverse-BFGS implementation and projected bounded variant.
   Consequently `optim.method` and arbitrary R `control` list semantics do not
   apply.

3. **Random-number stream**
   The supplied C++ uses separate `std::mt19937` engines for integer and uniform
   draws. The Fortran port uses the standard `random_number` generator seeded
   deterministically from `options%seed`. Identical seeds therefore do not
   imply identical cross-language population trajectories.

4. **Numerical derivatives**
   The public behavior is retained, but the Fortran implementation uses modern
   central finite differences rather than reproducing every legacy Gill/Murray
   derivative-estimation helper in `gradient.cpp` line for line.

5. **Parallel/R cluster execution**
   R's cluster, load-balancing, `share.type`, and environment-based parallel
   evaluation are orchestration features rather than the optimization kernel;
   they are omitted from the standalone FPM library.

6. **Project/output files**
   The obsolete/output-oriented `project.path`, `output.path`, and population
   dump machinery are omitted. Diagnostics are returned in `genoud_result`.

7. **`transform=TRUE` R wrapper**
   The R-specific transformed-return convention is not reproduced. The same
   calculations can be expressed directly in the Fortran objective callback.

8. **Custom R MemoryMatrix evaluator**
   The built-in duplicate cache is implemented. The R callback used to batch
   or externally parallelize MemoryMatrix evaluation is not reproduced.

The original package tree is included unchanged under `original/` so these
choices are auditable.
