# Porting notes

## Array and scalar behavior

The upstream R package vectorizes elementary error functions. In Fortran,
`se`, `ae`, `ape`, `sle`, and `ll` are `elemental`, so the same procedure can
be called with a scalar or any conformable array.

Aggregate functions require equal-length, nonempty vectors. R may recycle
vectors after a warning; the Fortran API rejects unequal shapes and returns
NaN with `metrics_invalid_size` when a status argument is supplied.

## Missing values

The upstream package does not provide an explicit missing-value policy. This
port does not skip NaNs. IEEE NaN or infinity values propagate through the
formulas, matching ordinary arithmetic behavior.

## Labels and retrieval lists

Classification error and accuracy have integer, double-precision real, and
character overloads. Information-retrieval set functions have the same three
overloads. Ragged `mapk` inputs use arrays of the small wrapper types
`integer_vector`, `real_vector`, and `string_vector`.

## Numerical edge cases

Zero denominators intentionally produce IEEE NaN or infinity where the R
formulas do. Invalid vector sizes and invalid binary labels are additionally
reported through optional status arguments.
