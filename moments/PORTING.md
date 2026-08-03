# Porting notes

## Data representation

R vectors are `real(dp) :: x(:)`. R matrices and data frames are represented by
`real(dp) :: x(:, :)`, with variables in columns. Allocatable function results
use Fortran's normal lower bound of one, so moment order `k` is stored at index
`k+1`.

## Missing values

IEEE NaNs are treated as missing values when `na_rm=.true.`. With
`na_rm=.false.`, a missing value produces an IEEE NaN result. The four test
routines remove nonfinite observations before applying their sample-size rules,
matching the practical intent of the upstream `complete.cases` calls while also
avoiding undefined arithmetic on infinities.

## Hypothesis-test alternatives

The original one-sided implementations associate R's `alternative="less"`
with the upper normal tail and `alternative="greater"` with the lower tail.
The Fortran port preserves that numerical convention for compatibility.

## Corrected cumulant recurrence

`all.cumulants` 0.14.1 initializes the first cumulant from the first central
moment, which is zero, instead of from the first raw moment, which is the mean.
That makes cumulants of order three and higher incorrect for nonzero-mean data.
The Fortran default implements the standard recurrence. Pass `legacy=.true.`
to reproduce the upstream values exactly.

## Error handling

Scalar descriptive functions return IEEE NaN for empty or degenerate input.
Hypothesis tests return explicit status values rather than raising an R error or
constructing an `htest` object.
