# neldermead-fortran

Modern Fortran/FPM translation of the computational algorithms in the R package
`neldermead` 1.0-13.

The package provides:

- `fminsearch`: MATLAB-compatible variable-shape Nelder-Mead defaults.
- `fminbnd`: Box's constrained complex method with bound and optional nonlinear
  positive inequality constraints.
- `neldermead_search`: configurable variable, fixed, or Box simplex optimization.
- `fmin_gridsearch`: direct grid evaluation.
- Native simplex construction/utility routines that replace the required
  computational subset of the R `optimbase`/`optimsimplex` dependencies.

## Build

```sh
fpm build
fpm test
```

GNU Fortran strict validation without FPM:

```sh
scripts/test_gfortran.sh
```

On Windows with gfortran:

```bat
scripts\test_gfortran.bat
```

See `API.md`, `TRANSLATION_COVERAGE.md`, and `VALIDATION.md`.
