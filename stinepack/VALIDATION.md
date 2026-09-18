# Validation notes

The maintained Fortran sources were compiled with GNU Fortran 14.2.0 using
Fortran 2018 mode, pedantic diagnostics, warnings, implicit-interface warnings,
line-truncation errors, and runtime checking.

The deterministic test program `test/test_stinepack.f90` covers:

- parabola slope recovery for a quadratic fixture;
- scaled and unscaled Stineman slope reference values;
- explicit-slope interpolation with a two-point fixture;
- scaled Stineman, Stineman, and parabola interpolation reference values;
- unique-prefix method matching and rejection of an ambiguous prefix;
- rejection of non-increasing knots and conflicting `yp`/`method` arguments;
- tiny endpoint extrapolation tolerance and NaN outside the allowed domain;
- vector gap filling, preservation/removal of unresolved edge NaNs;
- columnwise matrix gap filling and row omission.

An independent direct implementation of the formulas retained in `upstream/R`
was used as an additional deterministic comparison. With fixed random seeds,
120 estimated-slope interpolation cases (1,618 finite values) and 60
explicit-slope cases (674 finite values) matched the Fortran outputs exactly in
the validation environment (reported maximum absolute difference 0.0).

No test is evidence of complete R object/interface compatibility. The coverage
and documented differences in `README.md` and `API_COVERAGE.md` remain the scope
statement.

## Tool availability in the translation environment

The sandbox used for this translation did not contain the `fpm` or `fprettify`
executables, and its command-line network environment could not download them.
Accordingly, the literal `fpm build`, `fpm test`, and `fpm clean --all` commands
could not be executed here. The manifest was parsed as TOML and checked for
coverage consistency, while the exact library, test, and example source targets
were compiled and run directly with gfortran. Build products were removed
manually before packaging. Source formatting and line lengths were checked
without `fprettify`.
