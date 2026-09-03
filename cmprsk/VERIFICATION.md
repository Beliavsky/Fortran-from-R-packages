# Verification

Verification was performed on the maintained Fortran translation in this directory.

## Compiler regression

Compiler: GNU Fortran 14.2.0.

The library, six deterministic test programs, and the example were compiled directly with:

```text
-std=f2018 -Wall -Wextra -Werror=implicit-interface -fcheck=all -fbacktrace
```

Interface-compatible development stubs for the already-existing `rfortran-core` and `rfortran-linalg` APIs were kept outside the package tree. They were used only because those sibling directories are not mounted in this execution environment. The dependency names and called APIs were separately checked against the current target repository.

Result:

```text
test_api: PASS
test_crr_fit: PASS
test_crr_kernels: PASS
test_cuminc: PASS
test_gray_stratified: PASS
test_summary_timepoints: PASS
```

The example also compiled and ran successfully.

## Direct upstream-native parity

During development, the supplied original fixed-form routines were compiled outside the deliverable and evaluated on deterministic tied-time datasets. The modern translation was then evaluated on the identical inputs.

Exact-to-printed-precision comparisons were obtained for:

- `cinc`: all step times, cumulative-incidence estimates, and variances.
- `crstm`/`crst`: Gray scores and covariance, including a three-group/two-stratum case.
- `crrfsv`: Fine-Gray objective, score, and Hessian.
- `crrvv`: information and sandwich-meat matrices.
- `crrsr`: score residuals at each modeled-cause failure time.
- `crrfit`: baseline cumulative subdistribution-hazard jumps.

A complete Newton/Armijo fit driven by the original native routines produced, on the deterministic time-varying test problem:

```text
coef = [ 0.18003069023334850,
         1.6361297654789915,
        -1.2506868951378158 ]
log pseudo-likelihood = -8.0022561401616557
null log pseudo-likelihood = -8.6524231406763423
```

The high-level Fortran `fit_crr` reproduces these values along with the original robust covariance and baseline jumps; `test_crr_fit` also verifies the same result after shuffling input rows.

## Source-policy audit

`tools/audit_source.py` checks every maintained `.f90` file for:

- standard free-form line length;
- duplicate Fortran source content;
- executable semicolon-separated statements;
- `double precision`, `real*8`, `kind(0.0d0)`, and D-exponent literals;
- self-comparison NaN tests;
- separately declared dummy arguments;
- explicit `INTENT` or `VALUE` on every dummy argument;
- a trailing meaningful FORD `!!` comment on every dummy declaration;
- build products, binary libraries, modules, executables, caches, or ZIP archives inside the package tree.

Final audit: `AUDIT PASS: 16 Fortran files, 0 issues`.

## FPM availability

A real `fpm` executable is not installed in this execution environment. The required commands were explicitly attempted from the package root:

```text
fpm build
fpm test
fpm clean --all
```

Each returned exit status 127 (`fpm: command not found`). The exact transcript is retained in `FPM_ATTEMPTS.txt`. No fake or substitute `fpm` executable was used.

`fprettify` is also not installed in this environment, so it could not be run. The maintained source was nevertheless audited for free-form line length and formatting constraints and compiled in strict F2018 mode.
