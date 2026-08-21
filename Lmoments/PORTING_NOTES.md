# Porting notes

## Source mapping

- `R/Lmoments.R` + `src/Lmoments.cpp` -> `lmoments_core.f90`
- `R/QM.R` -> `lmoments_quantile_mixtures.f90`
- Hosking `SAMLMR`/`SAMLMU`, `SORT`, and `QUASTN` -> reference/utility logic in the modern Fortran implementation

## T1-L-moments

The Rcpp implementation builds recursive coefficient vectors that temporarily divide by zero at low order-statistic indices and then overwrites those entries. The Fortran port computes the same unbiased T1-L-moments directly from the standard order-statistic estimator

`E[X_(j:m)] = sum_i C(i-1,j-1) C(n-i,m-j) / C(n,m) * X_(i:n)`

and then applies the trimmed-L-moment alternating sum. This removes temporary infinities and makes the actual sample-size requirement explicit: estimating orders through `r` with trimming `(1,1)` requires at least `r+2` observations.

## Covariance

`lmom_cov` follows the current Rcpp `Lmomcov_calc` calculation. It evaluates the upper-triangular pair sum directly instead of materializing several `n x n` temporary matrices, reducing memory use from O(n^2) to O(n*rmax^2) while preserving the same arithmetic formula.

## Quantile mixtures

The package's CDFs use a fixed 49-step binary inversion (`j=-2,...,-50`). The Fortran implementation deliberately preserves that convention for parity. The quantile must be monotone for the returned density/CDF to define a valid distribution; like upstream, the port does not attempt to repair invalid coefficient vectors.

## R interface behavior omitted

R-specific `NA` removal, dimnames, list/class construction, warnings, and vector recycling are interface behavior rather than numerical algorithms and are not reproduced. Callers should remove missing values before passing arrays.
