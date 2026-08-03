# Testing

The test suite contains four programs:

- `test_moments`: direct moment formulas, optional tensor contractions, and
  magnitude adjustment.
- `test_sample_design`: `Q-MVSK`, `MM`, and `DC`, including simplex feasibility
  and objective improvement.
- `test_skew_t`: analytical skew-t moments, all six skew-t design methods, and
  fitted skew-t parameter estimation.
- `test_tilting_validation`: both tilting methods, tracking-error feasibility,
  delta/improvement consistency, and invalid-input status handling.

Validation commands:

```text
fpm test
```

or

```text
sh scripts/test_gfortran.sh
```

The shell script builds and runs strict runtime-checked and optimized
warning-as-error configurations, followed by all examples and the demo.
