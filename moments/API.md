# API

Import the complete public interface with `use moments`.

## Kinds and status values

- `dp = kind(1.0d0)`
- `MOMENTS_SUCCESS`
- `MOMENTS_INVALID_ARGUMENT`
- `MOMENTS_INSUFFICIENT_DATA`
- `MOMENTS_DEGENERATE_DATA`
- `MOMENTS_NONFINITE_DATA`

## Moment calculations

### `moment(x, order, central, absolute, na_rm)`

Returns a scalar for a vector or a vector of column results for a matrix.
Defaults match R: `order=1`, `central=.false.`, `absolute=.false.`, and
`na_rm=.false.`.

### `all_moments(x, order_max, central, absolute, na_rm)`

Returns every order from zero through `order_max`. Element `k+1` is order `k`.
The default maximum order is two.

### `skewness(x, na_rm)`

Population-moment skewness `mu3 / mu2**(3/2)`.

### `kurtosis(x, na_rm)`

Pearson kurtosis `mu4 / mu2**2`, not excess kurtosis.

### `geary(x, na_rm)`

Mean absolute deviation about the mean divided by the population standard
deviation.

## Moment transformations

### `raw2central(mu_raw)`

Converts a vector or column-wise matrix of raw moments to central moments.

### `central2raw(mu_central, eta)`

Converts central moments back to raw moments. `eta` is a scalar for a vector
and a vector of column means for a matrix.

### `all_cumulants(mu_raw, legacy)`

Calculates all cumulants through the maximum supplied order. The default uses
the standard raw-moment recurrence. Setting `legacy=.true.` reproduces the
upstream 0.14.1 initialization defect in which the first cumulant is set to the
first central moment instead of the mean.

## Hypothesis tests

The test functions return:

```fortran
type moments_test_result
   real(dp) :: statistic
   real(dp) :: transformed
   real(dp) :: p_value
   integer :: n
   integer :: alternative
   integer :: status
   character(len=48) :: method
   character(len=64) :: alternative_text
end type
```

Alternative constants are:

- `ALTERNATIVE_TWO_SIDED`
- `ALTERNATIVE_LESS`
- `ALTERNATIVE_GREATER`

### `agostino_test(x, alternative)`

D'Agostino transformed skewness test. Requires 8 through 46,340 complete
observations.

### `anscombe_test(x, alternative)`

Anscombe-Glynn transformed Pearson-kurtosis test.

### `bonett_test(x, alternative)`

Bonett-Seier test based on the mean absolute deviation and standard deviation.
As in R, `statistic` contains `tau`, and `transformed` contains the normal score.

### `jarque_test(x)`

Jarque-Bera statistic with the chi-square distribution on two degrees of
freedom. Its survival probability is evaluated exactly as `exp(-JB/2)`.
