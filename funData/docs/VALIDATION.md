# Validation

Validation was performed in the provided Linux sandbox with GNU Fortran 14.2.0.

## FPM command attempts

The requested literal FPM commands were attempted from the package root before packaging:

```text
$ fpm build
bash: line 1: fpm: command not found
exit=127

$ fpm test
bash: line 1: fpm: command not found
exit=127

$ fpm run --example basic_usage
bash: line 1: fpm: command not found
exit=127

$ fpm clean --all
bash: line 1: fpm: command not found
exit=127
```

The sandbox does not provide an `fpm` executable, so these attempts cannot be reported as successful FPM validation. The manifest was parsed as TOML and its source/test/example paths were audited directly.

## Direct GNU Fortran validation

To validate the same source/test/example units despite the missing FPM executable, two independent out-of-tree builds were run.

Runtime-checking build flags:

```text
-std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace
```

Optimized build flags:

```text
-std=f2018 -O2 -Wall -Wextra -Werror
```

For each configuration, the six library modules were compiled in dependency order, then both test programs and the example were compiled and executed. Results:

```text
All funData core tests passed
All funData simulation tests passed
integrals:    3.0000    4.0000
```

The simulation tests set a fixed Fortran random seed before stochastic operations. The deterministic core suite covers quadrature, 1D/2D/3D integration, regular and irregular norms, scalar products, regular/irregular conversion, missing-value interpolation, tensor products, orientation flipping including a regular-reference/irregular case, multivariate integration, and extraction behavior. The simulation suite covers eigenvalue/eigenfunction families, seeded univariate and multivariate simulation, sparsification, and Gaussian error addition.

## Release audits

Before creating the archive, automated checks verify:

- 20/20 coverage totals and one mapping entry per computational R export in the selected coverage basis;
- every mapped Fortran procedure is public and actually defined;
- all recorded `r_source` and `fortran_source` paths exist;
- one shared `dp = real64` definition and no legacy double-precision declarations or D-exponent literals;
- every dummy argument is declared separately with `INTENT` or `VALUE` and a trailing nonempty `!!` FORD comment;
- free-form ASCII Fortran with no maintained line longer than 132 characters;
- no disallowed code semicolons or self-comparison NaN tests;
- no system BLAS/LAPACK links, copied numerical dependencies, duplicate Fortran source files, compiler products, caches, or nested ZIP files.

Because `fpm clean --all` could not execute, compiler products from direct validation are kept outside the package tree and the tree is manually checked/cleaned before ZIP creation.
