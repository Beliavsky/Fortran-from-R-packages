# Translation notes

## Upstream

- Package: `rootSolve`
- Version: 1.8.2.4
- License: GPL >= 2
- Translation release: 0.1.0

## API mapping

| R routine | Fortran routine | Notes |
|---|---|---|
| `gradient` | `gradient` | Forward or centered numerical differences; same `max(abs(x)*pert,pert)` perturbation rule. |
| `hessian` | `hessian` | Scalar objective Hessian obtained by differencing numerical gradients. |
| `jacobian.full` | `jacobian_full` | Full numerical ODE Jacobian. |
| `jacobian.band` | `jacobian_band` | RootSolve-compatible compact band storage and grouped perturbations. |
| `multiroot` | `multiroot` | Newton-Raphson with scaled residual and step convergence tests. |
| `multiroot.1D` | `multiroot_1d` | Grouped banded finite-difference Newton iteration. |
| `uniroot.all` | `uniroot_all` | Grid sign-change scan followed by safeguarded Brent root refinement. |
| `stode` | `stode` | Dense/full or banded Newton steady-state solve; internal or user Jacobian. |
| `stodes` | `stodes` | CSR sparse Newton solve with discovered/user/grid sparsity. |
| `runsteady` | `runsteady` | ODEPACK DLSODE one-step integration with rootSolve's mean-absolute-derivative stopping criterion. |
| `steady` | `steady` | Dispatches to `stode`, `stodes`, or `runsteady`. |
| `steady.band` | `steady_band` | Banded `stode` convenience wrapper. |
| `steady.1D` | `steady_1d` | 1-D grid sparsity/banded/runsteady wrapper. |
| `steady.2D` | `steady_2d` | 2-D grid sparsity/runsteady wrapper. |
| `steady.3D` | `steady_3d` | 3-D grid sparsity/runsteady wrapper. |

## Numerical implementation details

### Finite differences

The upstream perturbation rule is preserved:

`delta_i = max(abs(x_i)*pert, pert)`.

`jacobian_band` and `multiroot_1d` perturb nonoverlapping band columns together, matching the optimization used by rootSolve's original banded code.

### Dense steady-state solver

The upstream `dSteady` algorithm is Newton-Raphson.  The translation preserves its two convergence tests:

1. `max(abs(f_i)/(rtol_i*abs(y_i)+atol_i)) <= 1`, and
2. maximum absolute Newton step <= `ctol`.

The original LINPACK `dgefa/dgesl` call is represented by a native pivoted Gaussian-elimination solve in the new rootSolve layer.

### Sparse steady-state solver

The original package contains Yale sparse direct solvers and SparseKit ILUT/ILUTP paths.  The translation preserves the sparse Jacobian target and grouped differencing but uses a modern CSR backend:

- caller-supplied `rowptr/colind`, automatically discovered sparsity, or 1-D/2-D/3-D grid sparsity;
- greedy column coloring for grouped finite differences;
- pivoted dense solve for `n <= 400`;
- diagonally preconditioned BiCGSTAB for larger systems.

The `spmethod` field is retained for API metadata, but this v0.1 implementation does not reproduce the instruction-for-instruction Yale/ILUT/ILUTP factorization kernels.  This is the principal algorithmic substitution in the translation.

### `runsteady`

This path deliberately reuses DLSODE from the free-format `deSolve` translation.  As in upstream rootSolve, DLSODE is called with `ITASK=2` (one internal step at a time).  After each step the RHS is reevaluated; convergence occurs when

`sum(abs(dy))/n < stol`.

If ODEPACK returns the excessive-accuracy code `-2`, absolute and relative tolerances are multiplied by ten and integration resumes, matching rootSolve's C driver.

### PDE wrappers

The 1-D/2-D/3-D sparsity maps use cell-major storage with species contiguous within each grid cell.  Each equation is connected to all species in the same cell and to the same species in immediate coordinate neighbors, with optional cyclic boundaries.  This is the structural pattern assumed by rootSolve's `sparse1d`, `sparse2d`, and `sparse3d` helpers.

## Deliberate R-interface omissions

The following are not numerical algorithms and are omitted:

- S3 `plot`, `image`, `subset`, and `summary` methods;
- R variable names/dimnames/class attributes;
- R `...`, list and environment handling;
- dynamic loading of user C/Fortran DLL callbacks;
- R forcing-table and output-variable plumbing.

Fortran callers provide procedure callbacks directly.

## Free-format status

All compiled Fortran in the release, including the vendored deSolve dependency, uses `.f90` free source form.  No fixed-form `.f`, `.for`, or `.f77` files are included.
