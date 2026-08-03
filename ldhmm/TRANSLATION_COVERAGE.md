# Translation coverage

## Included computational areas

- Complete `ecld` distribution calculations
- HMM model construction and parameter transformations
- Transition-matrix initialization and stationary probabilities
- Likelihood, forward/backward, local decoding, and Viterbi decoding
- Conditional distributions and pseudo-residuals
- State, density, and volatility forecasting
- Theoretical, decoded, and observation-conditioned statistics
- State and observation simulation
- Absolute-return ACF and moving-average utilities
- Maximum-likelihood optimization
- Generic price-to-log-return conversion

## Excluded areas

- `ldhmm.plot_spx_vix_obs`
- FRED downloading
- `xts`/`zoo` data classes and endpoint aggregation
- YAML configuration and package-data location helpers
- RData sample-object loading
- Package-specific CSV/archive readers
- Parallel R execution wrappers; Fortran routines are serial and reentrant
