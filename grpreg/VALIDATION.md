# Validation

Validation was performed on 2026-09-15 with GNU Fortran 14.2.0.

## FPM availability

The required FPM commands were explicitly attempted from the package root:

```text
$ fpm build
bash: line 1: fpm: command not found
exit status: 127

$ fpm test
bash: line 1: fpm: command not found
exit status: 127

$ fpm clean --all
bash: line 1: fpm: command not found
exit status: 127
```

FPM is not installed on this validation host, so these commands cannot honestly be
reported as successful. `fprettify` is also unavailable (exit status 127). The same
manifest source set, test, example, and independent parity program were therefore
compiled directly with gfortran in clean build directories. No compiler product was
written into the package source tree.

## Strict runtime-checked build

All maintained sources were compiled with:

```text
gfortran -std=f2018 -O0 -fcheck=all -fbacktrace \
  -Wall -Wextra -Werror -pedantic
```

Results:

```text
All grpreg deterministic tests passed.
Number of lambda values: 20
Smallest-lambda RSS:   7.2867E-14
```

The test suite exercises grouped lasso/MCP/SCAD, GEL, composite MCP, group bridge,
Gaussian/binomial/Poisson regression, Cox regression, coefficient/prediction/log-
likelihood/residual APIs, Gaussian/binomial/Poisson and Cox cross-validation, model
selection, mFDR, B-spline and natural-spline expansion, Cox survival/hazard/median
prediction, and the nonlinear data generator.

## Optimized build

The full test/example/parity sequence was repeated with:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -pedantic
```

All tests and the example passed. One `-O2` maybe-uninitialized diagnostic in the
mFDR working-weight allocation was resolved by making the allocation unconditional;
`-Werror` was retained.

## Independent numerical parity

`tools/parity_driver.f90` writes numerical results which are independently recomputed
by `tools/check_parity.py` with NumPy/SciPy rather than by another copy of the Fortran
algorithms. Both checked and optimized builds gave:

```text
gaussian_group_lasso: max_abs_diff=1.375e-11 limit=5.0e-09
binomial_near_mle: max_abs_diff=1.204e-07 limit=5.0e-06
poisson_near_mle: max_abs_diff=1.918e-08 limit=5.0e-06
cox_near_mle: max_abs_diff=2.774e-08 limit=5.0e-06
cox_deviance: max_abs_diff=0.000e+00 limit=5.0e-06
bspline: max_abs_diff=1.110e-16 limit=5.0e-12
Independent grpreg NumPy/SciPy parity checks passed.
```

The Gaussian reference includes a genuinely penalized grouped-lasso fit, not only a
near-unpenalized limit. Binomial, Poisson, and Cox checks independently optimize their
likelihoods with SciPy. The B-spline check uses SciPy's spline design matrix.

## Static and packaging audit

`python3 tools/audit_source.py` reports:

```text
SOURCE AUDIT PASSED: 15 Fortran files; 19/19 coverage mapping consistent.
```

The audit checks the single `dp` definition, line length, dummy `INTENT`/FORD comments,
legacy real-kind/D-exponent use, executable semicolon-separated statements,
self-comparison NaN idioms, duplicate maintained source files, build/archive products,
forbidden vendored dependency trees, and agreement of README/fpm.toml coverage counts.

A final release archive is additionally extracted into a fresh directory and the same
strict and optimized test/example/parity sequences are rerun against the extracted
source before release.
