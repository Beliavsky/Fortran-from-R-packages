# Validation

The maintained Fortran sources are validated independently of the retained upstream snapshot.

## Deterministic regression coverage

`test/test_earth.f90` checks:

- identity response contrasts;
- binomial-pair count expansion, row indexing, false/true ordering, and earth's zero-zero rule;
- weighted, multiple-response exact linear regression;
- coefficients, ordinary residuals, leverages, model-matrix/prediction consistency, deviance, and extractAIC/GCV semantics;
- exact piecewise hinge regression with an irrelevant predictor and `evimp` normalization;
- numerical parity with the upstream saved `trees` regression in `upstream/tests/test.earth.Rout.save` (selected basis, coefficients, RSS, R-squared, and GCV);
- exact degree-two hierarchical interaction construction.

`example/basic_earth.f90` fits and predicts a deterministic one-hinge regression.

## Compiler validation

The direct validation path uses GNU Fortran with strict Fortran 2018 diagnostics and runtime checking:

```text
gfortran -std=f2018 -pedantic -Wall -Wextra -Wimplicit-interface -Werror \
         -fcheck=all -fbacktrace -O0 -g
```

All maintained library sources are compiled first, followed by the test and example programs.

## FPM availability

The target commands remain:

```text
fpm build
fpm test
fpm run --example basic_earth
fpm clean --all
```

On the validation host used to create this archive, `fpm` is not installed and external DNS/package downloads are unavailable. Therefore those FPM commands cannot truthfully be reported as executed here. The manifest is dependency-free and the same source/test/example set is instead compiled and run directly with gfortran. The missing FPM executable is recorded as an environment limitation, not counted as a package test pass.

## Static policy audit

`validation/audit_source.py` checks maintained Fortran for free-form-only source, duplicate source content, disallowed semicolon statements, legacy real kinds/D exponents, self-comparison NaN idioms, the single `dp` definition, dummy-argument INTENT/VALUE plus FORD comments, build products, coverage consistency, and forbidden system-library linkage or shared compatibility files.
