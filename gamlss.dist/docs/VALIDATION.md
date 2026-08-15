# Validation

The v0.3.0 tree was rebuilt from an empty module directory with GNU Fortran
14.2.0 using:

```text
gfortran -std=f2018 -Werror=implicit-interface -fcheck=all -O0
```

No unlimited-free-line-length option is needed; every Fortran source, test, and
example line is within the standard 132-column free-form limit.

## Tests

Expected clean-build output:

```text
test_core: PASS
test_fit: PASS
test_v02: PASS
test_v03: PASS
```

`test_v03` adds checks for:

- independent scalar values for SN1, SN2, GT, ex-Gaussian and Pareto laws
- ST3C/ST3 numerical equivalence
- SST and SN1 CDF/quantile round trips
- exact double-binomial normalization
- PIG2/PIG parameter mapping
- zero-inflated/adjusted PIG, Sichel, BB, BNB and Zipf identities
- GAF/gamma and NBF/negative-binomial parameter mappings
- ZINBF zero-mass identity
- ZAPIG and ZASICHEL numerical normalization
- intercept-only Pareto-I recovery with `fit_gamlss`
- DBI recovery with `fit_dbi`
- zero-inflated beta-binomial recovery with `fit_zibb`

The v0.2 real-order modified-Bessel-K comparison against SciPy remains part of
the retained validation record; maximum absolute error in `log(K)` on that grid
was approximately `2.22e-9`.

## Examples

`basic`:

```text
BCT 95% quantile:   14.43630
NBI P(Y <= 4):       0.80679
```

`v02_extended`:

```text
ST3 90% quantile:           2.419381
SICHEL P(Y = 4):            0.081910
SHASH CDF at x = 0.7:       0.660794
```

`v03_remaining`:

```text
GT 90% quantile:                 1.204457
Double-binomial P(Y=3):         0.363325
Zero-adjusted PIG P(Y=0):       0.200000
Pareto-II 75% quantile:         1.816355
```
