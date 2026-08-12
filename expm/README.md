# expm-fortran

Modern Fortran/FPM translation of the computational code in the R package
**expm 1.0-0** (CRAN release dated 2024-08-19).

The library provides matrix exponentials, logarithms, square roots, integer
matrix powers, balancing, Frechet derivatives, exponential condition numbers,
and direct computation of `exp(t*A) v`.  R S3/S4/Matrix dispatch, printing,
plotting, localization, and `.Call`/`.Fortran` glue are intentionally omitted.

## Main API

Use the convenience module:

```fortran
use expm_module
```

The principal routines are:

- `expm(A [, balancing])` / `expm_higham08`: scaling and squaring with adaptive
  Pade orders 3, 5, 7, 9, and 13, plus optional LAPACK balancing.
- `expm_pade(A, order)`: translation of the R `R_Pade` scaling/Pade/squaring
  path.
- `expm_rbs(A [, degree, t])`: Roger Sidje/EXPOKIT-style diagonal Pade with
  scaling and squaring, corresponding to `PadeRBS`.
- `expm_taylor(A, order)`: scaled Taylor-series path corresponding to the
  package's legacy Taylor implementation.
- `expm_ward77(A [, order])`: trace shift, permutation/scaling balancing,
  Pade approximation, and reversal of the Ward-style preconditioning.
- `expm_almohy09(A [, p])`: translation of the package's configurable-degree
  scaling/Pade C kernel from `matexp_MH09.c`.
- `expm_eigen`, `logm_eigen`, and `expm_hybrid_eigen_ward`: spectral and hybrid
  counterparts of the R `R_Eigen`, `Eigen`, and hybrid paths.
- `matrix_power(A,k)`: nonnegative integer powers by binary exponentiation.
- `balance_real` / `balance_complex`: LAPACK `xGEBAL` wrappers returning the
  balanced matrix, scale vector, and active index range.
- `expm_frechet_sps`: Higham scaling/Pade/squaring Frechet derivative algorithm.
- `expm_frechet_block`: block-enlargement identity used by the original tests.
- `expm_cond_exact`, `expm_cond_1_est`, `expm_cond_f_est`: exact and estimated
  exponential condition numbers.
- `exp_at_v`: Sidje/EXPOKIT Krylov algorithm for `exp(t*A) v` without forming
  `exp(t*A)`.
- `sqrtm`: Schur-based principal matrix square root (complex result where needed).
- `logm`: inverse scaling and squaring with repeated Schur square roots and an
  atanh series.

Real and complex matrix exponentials are supported by the generic `expm`.
`sqrtm` and `logm` return complex matrices so negative-real spectral cases are
represented rather than silently discarded.

## Building with FPM

The manifest links BLAS and LAPACK:

```text
fpm build
fpm test
fpm run --example basic
```

A BLAS/LAPACK implementation visible to the linker is required.  The source is
standard free-form Fortran 2018 and does not require preprocessing.

## Validation

The release was compiled with GNU Fortran 14.2.0 using both an optimized build
and a warning-heavy bounds-checked build:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface -Werror
-std=f2018 -O2
```

Seven regression programs pass.  Highlights include:

- Ward (1977) 3x3 reference exponential.
- All translated exponential paths (`Higham08`, Pade, Taylor, RBS,
  Al-Mohy/Higham kernel, Ward).
- The expm package's 10x10 exact condition-number example:
  - 1-norm condition: `137.455837652872`
  - Frobenius condition: `566.582631819923`
- SPS Frechet derivative versus block enlargement.
- `exp_at_v` versus explicitly formed `exp(t*A) v`.
- Defective Jordan-block logarithm.
- Complex exponential and a square root with a negative eigenvalue.
- Deterministic matrix sizes 2 through 8 checking `exp(A)exp(-A)=I`,
  `log(exp(A))=A`, Frechet agreement, and square-root reconstruction.

FPM itself was not installed in the translation environment, so the FPM source
layout was validated with direct `gfortran` compilation and BLAS/LAPACK linking.

## Translation notes / intentional consolidation

The R package carries several historical implementations of essentially the
same scaling/Pade or Taylor computation (`Pade`, `PadeO`, `R_Pade`, etc.).
Rather than duplicate old and new loop bodies solely to preserve R method-name
aliases, this port exposes the distinct numerical kernels once through explicit
Fortran procedures.

The `sqrtm` implementation follows the same Schur-function principle as the R
code, using a complex Schur form and the triangular square-root recurrence.
`logm` retains inverse scaling-and-squaring and repeated square roots, but uses
an atanh-series final evaluation instead of reproducing the R
`logm.Higham08` rational-quadrature tables line-for-line.  This is the main
algorithmic consolidation in v0.1.0.

See `TRANSLATION_NOTES.md` for a source-to-module map and `LICENSES.md` for
licensing details.
