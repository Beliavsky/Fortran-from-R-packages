# Validation

The checked and optimized test suites cover:

1. BJ/GRS numerical equivalence, alpha-test output ranges, strong alternatives,
   and singular-design handling.
2. F2, Huberman-Kandel, and Kempf-Memmel frontier tests and finite-sample guards.
3. Reproducible Gungor-Luger sign-flip tests and rejection under a strong joint
   alternative.
4. Subseries Cauchy result naming, p-value ranges, B-draw reproducibility, and
   detection of an alpha alternative.
5. All twelve simulation presets, deterministic seeds, the exact GARCH(1,1)
   recursion, and invalid-input handling.
6. Normal, Student-t, and F probability-function fixtures.

`test_classical_alpha` additionally verifies the upstream documented identity
that Britten-Jones and GRS produce the same statistic and p-value.
