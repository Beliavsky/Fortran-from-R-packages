# Validation

The test suite checks:

1. Anderson-Darling statistics against direct evaluations of the upstream R formulas.
2. Marsaglia CDF and p-value fixtures in both source-compatible and clamped modes.
3. The uniform-data API and typed result fields.
4. The arbitrary-CDF callback API using a standard-normal transformation.
5. Empty, out-of-range, exact-endpoint, and clipped-endpoint behavior.

Reference fixtures were evaluated independently from the R formulas in
`R/ad.test.statistic.R` and `R/ad.test.pvalue.R` retained under `upstream/`.
