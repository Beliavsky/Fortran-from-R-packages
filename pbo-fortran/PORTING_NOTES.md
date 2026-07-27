# Porting notes

## Array orientation

As in the R package, input is `T x N`: observations are rows and candidate
strategies are columns. CSCV partitions contiguous blocks of rows and chooses
all `S choose S/2` half-sample combinations in lexicographic order.

## Callback interface

R accepts an arbitrary function returning one performance value per column.
Fortran uses an explicit procedure interface. The output vector is allocated by
the caller and must have one element per strategy. Internal procedures may be
used when a callback needs to capture parameters such as a risk-free rate.

## Ranking compatibility

The selected in-sample strategy is the first maximum, matching `which.max`.
The corresponding OOS score receives the average rank among exact ties,
matching the default `rank` behavior. Ranks run from 1 to N. When the normalized
rank is one, its positive-infinite logit is replaced by `inf_sub` exactly as in
the package.

## Upstream regression behavior

The R code constructs a two-column data frame named `Rn`, `Rbn`, then calls:

```text
lm(rn_pairs)
```

R's data-frame formula method interprets this as `Rn ~ Rbn`, not the degradation
plot's natural `Rbn ~ Rn` direction. The code then stores coefficient 1, the
intercept, in a field named `slope`, and coefficient 2, the slope, in a field
named `intercept`.

For compatibility, the Fortran fields `slope`, `intercept`, and `adjusted_r2`
reproduce that behavior and the original significant-digit rounding. The
additional fields `degradation_intercept`, `degradation_slope`, and
`degradation_r2` provide the directly interpretable regression:

```text
selected OOS performance = intercept + slope * selected IS performance
```

## Stochastic dominance data

The R dominance plot labels `CDF_all - CDF_selected` as `SD2`. That quantity is
a pointwise CDF difference rather than an integrated second-order stochastic
dominance measure. `dominance_curve` preserves it as `sd2_difference` and also
returns `integrated_difference`, calculated by cumulative trapezoidal
integration.

## Stronger validation

The R implementation relies on downstream errors for several invalid inputs.
The Fortran implementation reports unsuccessful results for odd subset counts,
subset counts that do not divide T, nonpositive `inf_sub`, and NaN metric
outputs. The threshold is a typed real value; the R package's accidental
acceptance of character thresholds is not reproduced.

## Numerical and presentation exclusions

No plotting was translated. No numerical content is lost: all plot inputs are
returned as arrays. The optional R parallel adapter is infrastructure rather
than a distinct algorithm and is not included.
