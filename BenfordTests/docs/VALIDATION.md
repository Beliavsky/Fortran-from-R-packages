# Validation

The Fortran translation is compiled and tested with GNU Fortran using:

```text
-std=f2018 -Werror=implicit-interface -Werror=trampolines -fcheck=all -O0
```

For the deterministic first-digit sample with counts
`(3,2,1,1,1,1,1,1,1)`, independent Python/SciPy calculations give:

| Statistic | Reference |
|---|---:|
| Pearson chi-square | 1.0958006038990171 |
| KS | 0.4005771787895563 |
| Chebyshev distance | 0.17677309040006986 |
| Euclidean distance | 0.30952309639746184 |
| Freedman-Watson U-square | 0.01579298003346328 |
| Mean-digit a-star | 0.11566974390945735 |
| JP-square | 0.9102630309457205 |
| chi-square asymptotic p-value | 0.9975676371883492 |

These are embedded as regression-test oracles. The joint-digit test with all
nonzero eigencomponents is also checked against the Pearson statistic, as
required by the multinomial covariance pseudoinverse identity.

The simulation tests are checked for support/shape rather than exact random
streams because Fortran does not standardize `random_number` sequences across
compilers.
