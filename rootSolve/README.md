# rootSolve-fortran

Modern free-format Fortran/FPM translation of the computational core of the R package `rootSolve` 1.8.2.4.

## Scope

The library provides native Fortran equivalents of the exported numerical routines:

- `gradient`
- `hessian`
- `jacobian.full`
- `jacobian.band`
- `multiroot`
- `multiroot.1D`
- `uniroot.all`
- `stode`
- `stodes`
- `runsteady`
- `steady`
- `steady.band`
- `steady.1D`
- `steady.2D`
- `steady.3D`

Plotting, image/subset/summary S3 methods, R dynamic-library loading, R expression/list plumbing, and graphics utilities are intentionally omitted.

## Design

The public API uses explicit Fortran procedure callbacks and derived result/option types.  Dense Newton methods use pivoted Gaussian elimination.  Banded finite-difference Jacobians preserve rootSolve's grouped perturbation strategy.  The sparse solver stores Jacobians in CSR form, colors nonconflicting columns so multiple variables can be perturbed in one function evaluation, and solves Newton systems with a pivoted dense solve for moderate systems or diagonally-preconditioned BiCGSTAB for larger systems.

`runsteady` uses ODEPACK DLSODE in one-step mode, exactly matching the key rootSolve algorithm: after each ODEPACK step it evaluates the model and stops when the mean absolute derivative is below `stol`.  The free-format deSolve translation is vendored as an FPM path dependency to supply the tested ODEPACK backend.

## Free-format requirement

Every Fortran file in this release is free-format `.f90`.  The original rootSolve fixed-form `.f` sources are deliberately not copied into the release tree; upstream metadata, R code, and C interface code are retained under `upstream/` for provenance.  The vendored deSolve backend is also entirely free-format `.f90`.

## Build

```text
fpm build
fpm test
```

The package has no C/C++ runtime dependency.  C files under `upstream/` are provenance only and are outside all FPM source directories.

## Minimal example

```fortran
program demo
  use rootsolve, only : dp, multiroot, root_result
  implicit none
  type(root_result) :: ans

  ans = multiroot(model, [1.0_dp, 1.0_dp])
  print *, ans%root
contains
  subroutine model(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = x(1)**2 + x(2)**2 - 1.0_dp
    f(2) = x(1)**2 - x(2)**2 + 0.5_dp
  end subroutine model
end program demo
```

The positive-start solution is approximately `(0.5, 0.866025403784...)`, matching the upstream `multiroot` example.
