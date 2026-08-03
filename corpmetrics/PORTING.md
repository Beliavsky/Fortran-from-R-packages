# Porting notes

## Scope

Every executable numerical function exported by corpmetrics 1.0 is translated:
`balsh`, `capm`, `ddm`, `fis`, `idm`, `insta`, and `loan`.

The R package has no native code, hidden numerical algorithms, plotting
functions, datasets, or nontrivial imported packages beyond basic `stats`
operations. R data frames, character formatting, and package documentation
objects are not computational algorithms and are not reproduced as runtime
objects.

## Preserved source behavior

### CAPM

Beta uses sample covariance divided by sample market variance. The risk-free
series is averaged independently, so its length need not equal the return
series, although the R documentation recommends equal lengths.

### Differential-growth DDM

The port intentionally uses the exact formula in `R/ddm.R`, including its
`DIV * (1 + G1)^(PER + 1)` terminal-dividend expression. It is not replaced by
a textbook alternative involving `G2` in the first terminal dividend.

### Fixed-income modified duration

For semiannual bonds, the original source divides Macaulay duration by
`1 + YTM`, not `1 + YTM/2`. The Fortran result preserves that formula.

### Investment analysis

NPV uses a potentially different cost of capital at each date. IRR ignores the
input cost vector after validating its length and solves for a single common
rate, exactly as the nested R function does.

### Loan tables

The recursive balance uses full precision. Interest, principal, and balance
columns are then rounded to two decimal places for output. The Fortran port
implements ties-to-even rounding through `round2_r`.

## Deliberate robustness improvements

- Zero-interest loans use the analytical installment `amount / periods`; the R
  expression produces `0/0`.
- IRR evaluation starts just above -100% to avoid division by zero at the
  original `uniroot` lower endpoint.
- Failed IRR brackets return `cm_root_not_bracketed` and IEEE NaN rather than
  terminating the process.
- Undefined ratios and singular valuation denominators return explicit status
  codes rather than relying on R warnings or infinities.
- Optional EPS and P/E values use logical availability flags instead of R `NA`.

These changes affect only invalid or limiting inputs, except where the source
formula is explicitly noted above and preserved.
