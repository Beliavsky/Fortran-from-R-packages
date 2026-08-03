# Porting notes

## R data frames and dates

The R function accepts a data frame and names of date/frequency columns. The
Fortran API accepts numeric arrays and integer period indices. This preserves
the model while avoiding dependence on R classes and supports irregular mixed
frequencies directly.

## Optimization

The upstream package combines `constrOptim`, `optim`, `maxLik`, and `numDeriv`.
The Fortran port uses self-contained BFGS and Nelder-Mead routines. Positivity,
stationarity, and beta-weight shape constraints are imposed by smooth
transformations rather than inequality constraints:

- `alpha >= 0`
- `beta >= 0`
- `gamma >= 0`
- `alpha + beta + gamma/2 < 1`
- weight-shape parameters `w1, w2 >= 1`

## Inference

The numerical Hessian and observation-level score matrix are calculated in
unconstrained coordinates. Covariances are then transformed to physical model
parameters with a numerical Jacobian. The port reports:

- inverse-Hessian covariance;
- sandwich covariance `H^-1 S'S H^-1`;
- OPG covariance with the same residual-kurtosis multiplier used upstream.

## Simulation corrections

Two clear edge cases in the original implementation are handled according to
the intended formulas:

1. When the low-frequency period contains multiple days, correlated covariate
   innovations use a standardized block return shock instead of R vector
   recycling.
2. In the positive-return branch of the realized-variance-dependent C++
   simulator, the upstream expression `r/r` is replaced by the intended
   squared standardized return `r^2/tau`.

These corrections are documented because exact reproduction of those coding
artifacts would not represent the stated model.

## Missing values

Unavailable leading MIDAS lags are represented by IEEE quiet NaNs. Likelihood
and variance-ratio calculations skip this leading unavailable region. Interior
missing values are not imputed.

## Omitted code

- plotting methods and weight plots;
- S3 print/predict dispatch (the numerical prediction is retained);
- R data objects and data-preparation scripts;
- R-specific validation of column names and `Date` classes.
