# Validation record

## Scope

This translation was validated on the maintained free-form Fortran sources in this
package. No source from `rfortran-core`, `rfortran-linalg`, `roll`, `rumidas`,
`rugarch`, or `maxLik` is copied into `dccmidas`.

The execution environment used for this translation does not contain FPM or local
checkouts of every sibling dependency. For direct compiler validation, small API
stand-ins were therefore compiled **outside** the package directory. They provide
only the documented procedure/type signatures needed to compile and exercise the
`dccmidas` sources and are not part of the archive. This direct validation checks
Fortran syntax, interfaces used by the package, array bounds/runtime behavior, and
the deterministic package-level numerical tests. It is not a substitute for the
final FPM integration test against the real sibling packages.

## Direct GNU Fortran validation

Compiler:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
```

The dependency stand-ins were compiled separately. The maintained package sources,
test, and example were then compiled with Fortran 2018, optimization, runtime
checking, warnings, and line truncation treated as an error. In abbreviated form:

```text
gfortran -c -std=f2018 -O2 -g -fcheck=all -Wall -Wextra \
  -Werror=line-truncation -Jmod -Imod src/*.f90

gfortran -std=f2018 -O2 -g -fcheck=all -Wall -Wextra \
  -Werror=line-truncation -Imod test/test_dccmidas.f90 ... -o test_dccmidas
./test_dccmidas
```

Result:

```text
All dccmidas tests passed.
```

The example was compiled with the same package-source flags and produced:

```text
Final conditional correlation:   0.833941
```

The deterministic test suite exercises determinant/inverse wrappers, RiskMetrics,
moving covariance, scalar and diagonal BEKK likelihood/matrix recursions, cDCC,
aDCC, DCC-MIDAS, asymmetric DCC-MIDAS, DECO, all five `cov_eval` loss choices,
QML standard errors, and high-level input validation.

## Static release audit

A source/package audit reported zero issues for:

- free-form `.f90` maintained source only;
- maximum source line length of 132 characters;
- no disallowed semicolon-separated executable statements;
- no self-comparison NaN tests;
- no `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals;
- every dummy argument declared by itself with explicit `INTENT` or `VALUE`;
- every dummy declaration carrying a meaningful trailing FORD `!!` comment;
- no duplicate maintained Fortran source files;
- no copied sibling dependency source trees;
- no object/module/executable/cache/ZIP build products in the package tree;
- all 22 `fpm.toml` function mappings resolving to existing upstream and Fortran
  source paths;
- README and manifest agreement at 22 of 22 mapped computational functions.

## FPM commands

The required FPM commands were explicitly attempted from the package directory:

```text
fpm build
fpm test
fpm run --example basic_dccmidas
fpm clean --all
```

All four attempts returned:

```text
fpm: command not found
```

Therefore this validation record does **not** claim that an FPM integration build
against the real sibling dependency checkouts was performed. After placing this
package at the repository root beside the dependencies named in `fpm.toml`, those
four commands should be rerun before treating the package as release-validated.
