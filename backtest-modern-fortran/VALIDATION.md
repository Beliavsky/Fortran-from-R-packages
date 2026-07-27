# Validation

Validation is performed by:

```text
make check
make release-check
```

The test suites cover:

1. R type-7 quantile values and exact bucket assignments.
2. Duplicate quantile-break rejection.
3. Global and period-by-period bucketing.
4. Group means, counts including NAs, and separate NA counts.
5. Strict trimmed-return endpoints.
6. Long/short turnover with exact hand-computed values.
7. Tri-bucket base weights and period normalization.
8. Two-period overlapping weights with exact hand-computed values.
9. Multi-period numeric backtests with exact low/high means and spreads.
10. Natural-portfolio turnover under full and partial constituent replacement.
11. Confidence-band finiteness and raw Sharpe ratios.
12. Total counts, marginal counts, cumulative bucket returns, and drawdowns.
13. Multiple signal and return columns using grouped inputs.
14. Numeric secondary-variable bucketing.
15. Logical universe filtering.
16. Demo, one-period CSV, overlapping CSV, and overlap example execution.

Exact compiler commands and final output are recorded in `TEST_RESULTS.txt` in
the packaged release.

No exact R object, print-layout, plot, or date/factor metadata equivalence is
claimed. The numerical definitions listed in `API_MAP.md` are the validated
translation target.
