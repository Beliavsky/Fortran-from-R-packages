# Validation

The release was compiled with gfortran 14.2.0 using:

    -std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all

and linked against BLAS/LAPACK.

All five permanent FPM tests pass:

1. `test_core`: independently known constrained optimum and objective.
2. `test_em`: EM does not increase the negative log-likelihood and preserves
   simplex normalization.
3. `test_preprocess`: row scaling invariance, log-likelihood input, all-zero
   column handling, and the one-positive-column trivial solution.
4. `test_svd`: full and low-rank paths agree on an exactly rank-3 positive
   matrix and the detected rank is 3.
5. `test_simulate`: generated dimensions, normalized likelihood range and RNG
   reproducibility.

The example program also compiles and runs under the same flags.

## Independent randomized optimization comparison

A separate development validation generated 120 random mixture problems with
6-24 rows and 2-7 columns, including cases with zero likelihood entries. Each
Fortran solution was compared with an independently formulated SciPy SLSQP
simplex-constrained minimization. All 120 Fortran runs converged.

- maximum absolute objective discrepancy: approximately `5.22e-14`
- maximum absolute mixture-coordinate discrepancy: approximately `3.70e-7`

The objective comparison is the stronger check in nearly non-identifiable
problems where several mixture vectors can have essentially identical
likelihood.
