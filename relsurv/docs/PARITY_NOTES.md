# v0.2.0 parity notes

## `rs.br`

The Fortran routine follows the upstream bridge construction: residuals are
scaled by event-specific variance, weights are proportional to `n.risk^rho`,
and the cumulative score process is converted to a Brownian bridge. Additive
model ties are collapsed using the upstream sum/sqrt(count) convention. Cox
ties retain individual residuals and multiply their weights by the tie
multiplicity before renormalization.

The maximum test uses the Kolmogorov bridge series. The Cramer-von Mises option
uses the upstream Watson-series probability calculation.

## `rs.zph` and `residuals.rsadd`

`rsadd_schoenfeld_residuals` translates the risk-set moments used by
`residuals.rsadd`, including the population and excess-hazard components and the
partial (`kvar`) covariance. Tied event times use averaged risk-set moments before
forming event residuals, as upstream does.

`rs_zph` implements identity, rank, log, and KM transforms and both upstream
variance modes: event-by-event inversion (`each`) and the summed covariance
scaling (`sum`).

## `epa`

`epa_smooth` implements the upstream adaptive segmentation of the event grid.
Interior evaluation uses the ordinary Epanechnikov kernel; left/right boundary
evaluation uses the asymmetric boundary polynomial; predictable left smoothing
uses the one-sided kernel. Native callers can evaluate the same kernel at an
arbitrary supplied target grid.

## `rsadd` EM and GLM branches

The EM implementation preserves cause codes 0/1/2, posterior unknown-cause
weights, weighted counting-process Cox M-steps, unsmoothed or left-Epanechnikov
smoothed baseline excess hazards, left truncation, and the missing-information
Fisher correction. Automatic bandwidth selection compares the p=0 EM baseline
against the translated Ederer-II target, matching the upstream algorithmic
criterion.

The grouped binomial and Poisson branches implement the custom population-
survival links/offsets rather than replacing them with generic GLMs.

## `years`

The Greenwood area variance is the upstream cumulative area formula, not a
delta-method approximation. YL2013 and YL2017 use their distinct left/right
integration conventions. When bootstrap replicate curves are supplied, column
sample variances are calculated for probabilities, integrated areas, and final
years estimates exactly as in the numerical aggregation stage of the R code.

## HLD/HMD import

The HLD parser filters the standard life-table rows, carries terminal qx values
where required, and converts annual death probabilities to daily hazards using
`-log(1-qx)/365.241`. Multiple HLD inputs can populate the race dimension. The
HMD parser enforces the 0:110 age grid and applies the upstream terminal-qx
handling before the same daily-hazard conversion.
