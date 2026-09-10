# Validation

The maintained Fortran implementation is self-contained and uses no external BLAS/LAPACK, R, Rcpp, or C++ numerical library.

## Compiler validation performed

The translation is validated with GNU Fortran 14.2.0 using both:

- strict debug checks: `-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -Werror=implicit-interface -fcheck=all -g`;
- optimized checks: `-std=f2018 -O2 -Wall -Wextra -Werror`.

The deterministic suite exercises:

- honest regression forests, missing-value routing, normalized forest kernels, tree/leaf/split/importance/merge utilities;
- grouped little-bag sampling, complete-group tree rounding, and debiased regression/causal/IV/causal-survival/local-linear variances;
- automatic regression hyperparameter tuning and resolved-option storage;
- multi-response, class-probability, and conditional-quantile forests;
- OOB nuisance forests for causal/IV/local-coefficient families, conditional treatment-variance scores, and binary-IV compliance forests;
- survival and censoring-adjusted causal-survival RMST/survival-probability workflows;
- ATE, calibration, BLP, RATE and orthogonal-score helpers;
- local-coefficient and multi-arm forests;
- local-linear OOB ridge-path selection, covariance-weighted ridge prediction, and local-linear residual split training;
- fixed-stage and automatically stopped boosted regression with OOB debiased-error diagnostics;
- all twelve deterministic causal benchmark generators and all six causal-survival generators, including scaling, correlation, supplied-predictor, Monte Carlo, and validation controls.

Both programs under `example/` are also built and run during release verification.

`tools/check_style.py` checks the maintained Fortran files for the requested dummy-argument documentation/INTENT rules and prohibited legacy real/semicolon constructs.

## FPM and formatting-tool availability

The build image does not contain an `fpm` executable. A prior attempt to retrieve a current standalone FPM release asset was blocked by the environment, so literal `fpm build` and `fpm test` cannot be claimed here. The included `fpm.toml` is parsed with Python's TOML parser, uses standard source/test/example discovery, has no dependencies, and the identical source/test/example graph is compiled directly with gfortran in both strict and optimized configurations.

`fprettify` is not installed in the build image. The source is therefore audited directly rather than claiming that `fprettify` ran.

## Upstream integrity

Before packaging the original supplied tree, all 242 entries in the upstream CRAN `MD5` manifest were verified. The generated upstream artifact `build/vignette.rds` is intentionally omitted from the deliverable because the translation ZIP must not contain build products. All retained upstream files are unchanged.
