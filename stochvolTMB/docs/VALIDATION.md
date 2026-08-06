# Validation

Six deterministic test programs cover:

1. simulation for all five models and standardized skew-normal moments;
2. finite Laplace objectives, positive latent Hessians, and latent-mode signal recovery;
3. Gaussian parameter optimization and objective improvement;
4. short fitting paths for Student-t, skew-normal, leverage, and skew-leverage models;
5. parameter simulation, predictive draws, and prediction summaries;
6. normal/Student-t/skew-normal probability functions and PIT residual moments.

The source tree is compiled in two configurations:

- checked: Fortran 2018, all warnings as errors, bounds checks, backtraces, and floating-point traps;
- optimized: `-O3 -march=native`, all warnings as errors.

The final archive is extracted into a separate directory and rebuilt in both
configurations before release.
