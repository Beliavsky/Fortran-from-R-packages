# Validation

Validation was performed with GNU Fortran 14.2 using Fortran 2018 mode,
runtime bounds/array checks, backtraces, and implicit-interface errors enabled.

The test suite covers:

- multinomial coefficients;
- reduction to an ordinary multinomial when every `theta=1`;
- a hand-solvable two-category interaction case (`dMM([1,1]) = 2/3`);
- independently computed 3-category normalizer and density references;
- weighted log likelihood and sufficient statistics;
- equality of the sufficient-statistic support likelihood plus the omitted
  multinomial constants with the full likelihood;
- independently computed expected sufficient statistics;
- Lindsey's Poisson device on a saturated two-category example with exact
  target `p=(0.4,0.6)` and `theta=1`;
- smooth/BFGS and Nelder-Mead likelihood refinement;
- differing-row-sum likelihood identity;
- multinomial and MB `gunter` support aggregation;
- bivariate `Lindsey_MB` coefficient construction/fit;
- MCMC output dimensions, non-negativity, and preservation of total count.

Independent numerical reference values used by `test_core` include:

```text
p = (0.2, 0.3, 0.5)
theta12 = 1.4
theta13 = 0.8
theta23 = 1.7

NormC(Y=3) = 1.853920000000001
P(Y=(1,1,1)) = 0.18486234573228613
weighted log likelihood = -13.16014610575452
```

These values were generated independently from a direct enumeration of the
multiplicative-multinomial formula, rather than from the translated routines.
