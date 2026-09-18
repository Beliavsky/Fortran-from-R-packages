# Validation

## Compiler checks performed

The maintained Fortran sources, deterministic tests, and example were compiled with GNU Fortran 14.2.0 using:

```text
-std=f2018 -pedantic -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -ffree-line-length-132
```

The real translated `stinepack` sources from the preceding sibling package were used during validation.

The sandbox had no `fpm` executable and no working command-line network path from which to install or clone it/dependencies. The required commands were attempted explicitly:

```text
> fpm build
bash: fpm: command not found
> fpm test
bash: fpm: command not found
> fpm clean --all
bash: fpm: command not found
```

Thus the literal FPM build/test/clean could not run in this environment. For direct compiler validation of the `forecast` calls, a temporary interface-compatible forecast shim was compiled outside the package. It supplied only the public types/signatures used by imputeTS and simple test implementations. The shim is not included in the package ZIP. Public signatures were checked against the current sibling `forecast` source in Beliavsky/Fortran-from-R-packages.

## Tests performed

`test/test_imputets.f90` passes and covers:

- `na_replace` / `na_remove`
- arithmetic, median, mode, geometric, and harmonic `na_mean` modes
- LOCF and NOCB edge policies
- simple, linear, and exponential moving-average weights
- linear, natural-cubic, and Stineman interpolation
- structural and ARIMA Kalman branches
- seasonal splitting and seasonal decomposition orchestration, including the upstream early-return fallback behavior
- random imputation bounds with a fixed processor RNG seed
- gap statistics including tie behavior
- `maxgap` restoration

The example also compiles and runs. The final direct gfortran regression run completed with `All imputeTS tests passed.`

## Source-policy audit

Before packaging, the maintained source is checked for: one real kind, forbidden D/d exponent forms, `double precision`/`real*8`, semicolon-separated statements, self-comparison NaN tests, missing dummy `INTENT`/`VALUE`, multi-dummy declarations, missing trailing FORD dummy comments, lines over 132 characters, duplicate Fortran sources, copied dependency trees, and build products.

`fprettify` is also unavailable in this sandbox (`fprettify: command not found`). Maintained free-form source was therefore checked mechanically for the requested line-length, declaration, operator-spacing, and related style constraints instead.

Because the actual `forecast` sibling could not be materialized in this sandbox, a final literal `fpm build && fpm test` should be run after extracting the package next to `forecast` and `stinepack` in the target repository.
