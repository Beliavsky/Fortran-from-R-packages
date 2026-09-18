# Validation

This translation was validated in the available Linux build environment with GNU Fortran and an independent Python/SciPy reference implementation. The package is intended to be FPM-compatible on Windows with gfortran as described in `README.md`.

## Toolchain

- GNU Fortran (Debian 14.2.0-19) 14.2.0
- Python 3.13.5
- NumPy 2.3.5
- SciPy 1.17.0
- `fprettify`: not installed in the validation environment
- `fpm`: not installed in the validation environment

## Required FPM commands

The requested FPM commands were attempted from the package root:

```text
fpm build          -> exit 127: fpm: command not found
fpm test           -> exit 127: fpm: command not found
fpm clean --all    -> exit 127: fpm: command not found
```

Because FPM was unavailable, these commands cannot truthfully be reported as successful. The same manifest-declared library sources, test program, and example were instead compiled directly with gfortran in clean external build directories. No FPM or compiler build products are included in the package archive.

## Strict debug/runtime-check build

The maintained sources were compiled with:

```text
-std=f2018 -O0 -g -fcheck=all -fbacktrace -Wall -Wextra -Werror -pedantic
```

Results:

```text
All poLCA tests passed.

log likelihood:   -1487.4508
class shares:      0.42976   0.57024
class 1 P(Y=2):   0.11271   0.17915   0.06904
class 2 P(Y=2):   0.79188   0.71988   0.77875
```

## Optimized build

The same tests, example, and parity driver were compiled independently with:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -pedantic
```

All deterministic tests passed and the example produced the same displayed values as the runtime-checked build.

## Independent numerical parity

`tools/parity_check.py` fits the same synthetic data by directly optimizing the observed-data likelihood with SciPy rather than reusing the translated EM algorithm. It also independently constructs the poLCA outer-product covariance and response-probability/class-share delta-method standard errors.

Both debug and optimized Fortran builds passed with the following maximum absolute errors:

```text
no-cov loglik                         5.09e-11
no-cov shares                         2.33e-08
no-cov response probabilities class 1 3.22e-08
no-cov response probabilities class 2 2.09e-08
no-cov class-share SE                 3.47e-10
no-cov response SE class 1            7.10e-10
no-cov response SE class 2            2.08e-10
regression loglik                      4.55e-12
regression beta                        2.11e-07
regression response probs class 1      4.50e-08
regression response probs class 2      7.68e-08
regression beta SE                     1.03e-08
regression class-share SE              1.05e-09
```

The parity script accepts the parity-driver output filename as its first argument, for example:

```text
python tools/parity_check.py path/to/fortran_parity.txt
```

## Source-policy audit

`python tools/source_audit.py` checks maintained Fortran source for:

- free-form line length at most 132 characters;
- disallowed semicolon-separated statements;
- legacy `double precision`, `real*8`, and D-exponent literals;
- self-comparison NaN tests;
- disallowed fast/finite-math options in maintained source;
- one dummy argument per declaration;
- explicit `INTENT` or `VALUE` and trailing FORD `!!` documentation for every dummy argument;
- duplicate maintained Fortran source content; and
- consistency of the 10/10 coverage metadata with `README.md`.

Final result:

```text
SOURCE AUDIT PASSED: 9 Fortran files; 10/10 coverage mapping consistent.
```

## Clean-tree checks

Before archiving, package-local compiler/build artifacts, caches, and temporary validation directories were removed manually because `fpm clean --all` could not run. The final archive is checked to contain one top-level `poLCA/` directory, no object/module/executable/cache/ZIP products, no copied dependency source tree, and no duplicate maintained Fortran sources.
