# Translation notes

## Source package

- R package: `expm`
- Source version: `1.0-0`
- Date: 2024-08-19
- Package license: GPL >= 2
- `src/matexp_MH09.c`: GPL >= 3

The original sources used for the translation are retained under `original/`.

## Source map

| Original code | Fortran translation |
|---|---|
| `R/expm.R`, `R/expm2.R` | `src/expm_matrix_functions.f90` |
| `src/matexp.f` | `expm_rbs` in `src/expm_matrix_functions.f90` |
| `src/matrexp.f`, `src/matrexpO.f`, `src/mexp-common.f` | `expm_taylor`, `expm_pade`, matrix helpers |
| `src/matexp_MH09.c` | `expm_almohy09` |
| `R/balance.R`, `src/R_dgebal.c` | `src/expm_linalg.f90` (`balance_real`, `balance_complex`) |
| `R/matpow.R`, `src/matpow.c` | `matrix_power` |
| `R/expm_vec.R` | `src/expm_action.f90` |
| `R/expmCond-all.R` | `src/expm_frechet.f90`, `src/expm_condition.f90` |
| `R/sqrtm.R` | `src/expm_logsqrt.f90` Schur square root |
| `R/logm.Higham08.R`, `R/logm.R` | `src/expm_logsqrt.f90` inverse scaling/squaring log |
| `src/expm-eigen.c`, `src/logm-eigen.c`, R spectral paths | `src/expm_eigen_methods.f90` |
| R `.Call`, S3/S4/Matrix dispatch | omitted |

## Numerical design

Fortran arrays and BLAS/LAPACK are column-major, so most matrix formulas map
naturally from R and the package's C/Fortran kernels.  All dynamically sized
work arrays use allocatable Fortran arrays and all LAPACK calls have explicit
interfaces.

### Matrix exponential

`expm_higham08` implements the package's Higham scaling/Pade/squaring structure
with the same low-order threshold table and the degree-13 coefficients used by
`expm.Higham08`.  Optional balancing performs the package's separate permutation
and scaling steps and reverses them explicitly.

The historical method family is exposed as distinct computational kernels,
without recreating R's method-string dispatcher.

### Frechet derivative and condition numbers

`expm_frechet_sps` is a direct translation of the R implementation of Higham's
2008 scaling/Pade/squaring Frechet algorithm.  The exact condition number builds
the full Kronecker representation exactly as the R code does.  The two
estimators follow the original 1-norm estimator and Frobenius power iteration.

### Matrix square root and logarithm

The R square-root code works from a real Schur decomposition and handles 1x1
and 2x2 blocks.  The Fortran version asks LAPACK for a complex Schur form; this
turns the quasi-triangular block recurrence into the standard upper-triangular
recurrence and naturally supports complex principal roots.

The R `logm.Higham08` uses inverse scaling/squaring followed by a rational
quadrature.  The Fortran version preserves inverse scaling/squaring and Schur
square roots but evaluates the final near-identity logarithm with

`log(X) = 2 * atanh((X-I)(X+I)^(-1))`

and the odd-power atanh series.  This avoids embedding a large table of
R-specific quadrature constants while retaining stable near-identity behavior,
including defective Jordan blocks.

## Not translated

The following are deliberately outside the computational Fortran library:

- R S3/S4 classes and `Matrix` coercion/dispatch.
- `.Call`, `.Fortran`, registration, and R memory-management code.
- printing, plotting, demos, gettext/localization, and package metadata logic.
- parallel R orchestration (none of the core dense numerical kernels require it).
- duplicate legacy loop bodies when they implement the same Pade/Taylor kernel.

## Dependencies

Only BLAS and LAPACK are required.  No R runtime or C shim is required.
