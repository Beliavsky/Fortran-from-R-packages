# Porting notes

## Scope

The attached roptim source contains the C++ wrapper and simulated-annealing
code, while Nelder-Mead, BFGS, CG, and L-BFGS-B are called from R's native
optimization library. A self-contained Fortran package therefore requires
numerical implementations in addition to translating the wrapper.

## Implementations

- Nelder-Mead uses reflection, expansion, contraction, and shrink operations
  controlled by `alpha`, `gamma`, and `beta`.
- BFGS stores a dense inverse-Hessian approximation and uses a safeguarded
  Armijo/Wolfe-style line search.
- CG supports three update rules corresponding to the package's `type=1..3`
  control and periodically restarts.
- L-BFGS-B uses the version 3.0 reverse-communication numerical kernel.
- SANN follows the package's logarithmic temperature schedule and Gaussian
  Markov proposal, with an explicit custom proposal callback.

## Differences from R/C++

- Arrays use normal Fortran one-based indexing.
- Bounds are passed directly to `roptim_minimize`.
- A custom SANN proposal receives `(current, candidate, scale, user_data)`;
  it is not overloaded onto a gradient method.
- Result messages are descriptive Fortran strings.
- Evaluation counts can differ from R because finite-difference and line-search
  implementations are native rather than calls into R's compiled library.
- Rcpp, Armadillo, S3 objects, R printing, and R random-number state are omitted.
- The Fortran intrinsic RNG is used for SANN; `control%seed` gives reproducible
  runs within a compiler/runtime but is not intended to reproduce R's RNG stream.

## GNU Fortran callback-interface compatibility

Version 0.1.1 routes user callbacks through module-level trampoline procedures
whose dummy procedures use the package's abstract interfaces explicitly. This
avoids a GNU Fortran diagnostic seen in some FPM/compiler combinations where a
host-associated optional procedure dummy was incorrectly reported as having an
implicit interface. The objective, gradient, proposal, and monitor callback
semantics are unchanged.
