# Porting Notes

## Scope mapping

The translation retains the computational content of the R and native C/C++
source:

- all 13 conditional-duration recursions
- all ten duration distributions present in the native likelihood code
- maximum-likelihood fitting and parameter inference
- score and covariance calculations
- simulation and forecasting
- transaction-to-duration conversion
- all four diurnal-adjustment methods
- residual transformations and numerical diagnostics
- all three model-specification tests

The six graphics-oriented exports are not rendered as plots. Their useful
numerical data are returned through summary, ACF, density, QQ, hazard, rolling,
and likelihood-profile result types.

## Solver substitutions

The original R package delegates optimization to `optim`, `optimx`, or
`Rsolnp`, and contains thousands of lines of generated analytical-score C/C++
code. The Fortran implementation uses:

- a self-contained bounded Nelder-Mead optimizer with restarts
- central or bounded finite-difference observation scores
- a numerical likelihood Hessian
- dense self-contained linear algebra for covariance calculations

This preserves the fitted objective and model equations, but iteration counts and
reported last digits can differ from an R optimizer run.

## Source inconsistencies resolved

The original package contains several discrepancies between R wrappers,
likelihood code, and simulation code. The Fortran port applies one documented,
internally consistent definition in each case.

1. **AMACD residual lags**
   The likelihood repeatedly used the first residual lag for every MA term,
   whereas the simulator indexed lag `j`. The Fortran recursion uses `i-j` for
   each residual lag.

2. **ABACD/AACD shape parameters**
   Native simulators and the likelihood disagree on whether asymmetry and power
   parameters are supplied once or once per lag. The R parameter-count and
   likelihood convention—one common `c`, `v`, `d1`, and `d2`—is used.

3. **AACD/ABACD defaults**
   The default vector is made consistent with the R start-value routine:
   `c=0.8`, `v=0.1`, `d1=1.1`, and `d2=1.1`.

4. **TACD/TAMACD layout**
   The likelihood's regime-major parameter layout is used consistently by
   filtering and simulation.

5. **SNIACD/LSNIACD indexing**
   The R parameter-count helper uses `min(0,p-1)`, while the start-value and
   intended recursion require `max(0,p-1)`. The intended positive lag count is
   used. Spline coefficients are shared by the first lag and optional scale
   multipliers are used for later lags.

6. **LSNIACD initialization**
   A native simulator initializes a log conditional mean with `exp(startMu)`.
   The Fortran code uses the mathematically consistent `log(startMu)` state.

7. **Generalized-F forced mean**
   The R density wrapper has a sign inconsistency in its forced-scale power.
   The likelihood/helper expression using `eta^(-1/gamma)` is used.

8. **q-Weibull formulas**
   The q-Weibull quantile wrapper refers to an undefined variable in one branch,
   and the documented defaults include both `q<1` and `q>1` cases. The Fortran
   implementation provides valid finite-support and heavy-tail formulas for both
   ranges.

9. **Mixture CDF wrappers**
   Duplicate/partial argument names in the R mixture CDF calls are replaced by
   explicit parameter positions.

## Initial conditions

As in ACDm, filtering begins each series or day segment with the sample mean for
at least the maximum lag. Simulation accepts explicit initial durations and
conditional means; otherwise it starts from unit conditional duration and uses a
burn-in.

## Random numbers

The Fortran library uses an explicit xorshift-family generator and its own normal
and gamma samplers. Seeds are reproducible within this project, but streams do
not match R's RNG exactly.

## Time representation

R POSIXlt values are replaced by integer calendar arrays and real seconds since
midnight. Diurnal routines accept numeric clock coordinates and optional integer
group IDs, making weekday or per-date aggregation explicit and portable.
