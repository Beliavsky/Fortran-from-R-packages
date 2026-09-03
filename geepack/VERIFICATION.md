# Verification

## Strict compiler regression

The maintained library was compiled from a fresh external build directory with
GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -fbacktrace
```

Because the target repository sibling packages are not mounted in this runtime,
small interface-compatible test stubs for `r_kinds`, `r_distributions`, and
`r_linalg` were kept **outside** the deliverable under `/mnt/data/geepack_build`.
The public dependency procedure signatures were verified against the current
repository before this build. No stub or dependency source is present in the
`geepack` directory or ZIP.

Result:

- 8 deterministic test programs passed.
- 1 example program compiled and ran successfully.
- Runtime array/bounds checking remained enabled for all tests.
- Implicit external procedure calls were treated as errors.

The deterministic coverage includes Gaussian and binomial GEE, mean and scale
models, all six ordinary working-correlation modes, custom `corp` AR(1)
coordinates, AJS/J1S/FIJ covariance paths, ordinal GEE with and without
association estimation, correlation/design helpers, QIC/CIC/QICC, COPY
relative-risk fitting, and Wald inference.

## Direct upstream parity checks

During translation, deterministic reference cases were compared with the
original package algorithms for link/variance/correlation construction and the
estimating-equation formulas. The regression tests retain independent numerical
references for OLS/logistic coefficients, QIC quantities, correlation-design
ordering, and finite-difference ordinal odds derivatives.

The 19 upstream files retained under `upstream/` were checked byte-for-byte
against the supplied `geepack-master.zip` extraction before packaging.

## Static/package audit

Final audit result: **20 Fortran files, 0 policy issues**.

Checks included:

- standard free-form line length (maximum 132 columns);
- every dummy argument separately declared with `INTENT` or `VALUE` and a
  meaningful trailing FORD `!!` comment;
- no executable semicolon-separated statements;
- no `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals;
- no self-comparison NaN idioms;
- no duplicate Fortran source contents;
- no build products or ZIP files inside the package;
- no copied `rfortran-core`, `rfortran-linalg`, BLAS, LAPACK, or ARPACK source;
- `fpm.toml` parses and contains only sibling path dependencies, with no system
  library links.

`fprettify` is not installed in the execution environment; source was therefore
kept manually in compatible free-form formatting.

## FPM command attempts

A real `fpm` executable is not installed in this runtime. The required commands
were nevertheless invoked from the package root immediately before packaging:

```text
fpm build
fpm test
fpm clean --all
```

All three returned shell exit status 127 (`fpm: command not found`). An attempt
to reach the upstream FPM release site from the container also failed because
outbound DNS is disabled (`Could not resolve host: github.com`). No fake or
substitute FPM executable was used. The exact command transcript is retained in
`FPM_ATTEMPTS.txt`.
