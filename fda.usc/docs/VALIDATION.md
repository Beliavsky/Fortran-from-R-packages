# Validation

Release validation uses GNU Fortran 14.2.0 in two direct-build modes:

1. runtime checking: `-std=f2018 -O0 -fcheck=all -fbacktrace`;
2. optimized: `-std=f2018 -O2`.

In each mode, all four deterministic test executables and the `basic` example are compiled and run. The four tests are `unit`, `statistics`, `extended`, and `generalized`. They exercise both the original translated kernels and the continued-translation depth/basis/regression/classification/test families. Tests do not require RNG equality with R.

The preparation environment does not provide an `fpm` executable. Literal `fpm build`, `fpm test`, `fpm run --example basic`, and `fpm clean --all` commands were attempted before packaging and returned exit 127 (`fpm: command not found`). These attempts are **not** reported as successful FPM runs.

To validate the same source graph despite that environment limitation, the maintained package sources are compiled directly with GNU Fortran against narrow external validation stubs matching the public interfaces used from `rfortran-core` and `rfortran-linalg`. The stubs are not part of the package or ZIP.

## Release results

Both direct-build modes passed:

```text
checked unit        -> All fda.usc unit tests passed
checked statistics  -> All fda.usc statistics tests passed
checked extended    -> All fda.usc extended tests passed
checked generalized -> All fda.usc generalized tests passed
checked example     -> L2 distance between curves 1 and 3: 1.154701

optimized unit        -> All fda.usc unit tests passed
optimized statistics  -> All fda.usc statistics tests passed
optimized extended    -> All fda.usc extended tests passed
optimized generalized -> All fda.usc generalized tests passed
optimized example     -> L2 distance between curves 1 and 3: 1.154701
```

Literal FPM attempts:

```text
fpm build                 -> exit 127: fpm: command not found
fpm test                  -> exit 127: fpm: command not found
fpm run --example basic   -> exit 127: fpm: command not found
fpm clean --all           -> exit 127: fpm: command not found
```

## Release audits

The final release process checks that:

- `fpm.toml` parses as TOML and records 179 of 179 mapped exported computational R functions;
- every mapping points to an existing public Fortran procedure;
- every retained `r_source` path exists under `upstream/R`;
- every Fortran dummy data argument has explicit `INTENT` or `VALUE`, is declared separately, and has a trailing nonempty FORD `!!` comment;
- maintained Fortran text is ASCII free-form and no line exceeds 132 characters;
- no disallowed code semicolons, self-comparison NaN idioms, `double precision`, `real*8`, `d0` literals, or package-local alternative real kinds are used;
- no system BLAS/LAPACK links are present;
- no dependency source is copied into the package;
- no duplicate Fortran source files are present;
- no object/module/executable/cache/ZIP artifacts are packaged;
- the final archive contains exactly one top-level `fda.usc/` directory.

Because FPM is absent, `fpm clean --all` cannot execute; the source tree is instead independently audited to confirm that it contains no compiler or FPM build output before archiving.
