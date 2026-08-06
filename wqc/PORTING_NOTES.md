# Porting notes

## Main translation choices

1. The R `lapply`/`sapply` result assembly is represented by derived types and
   allocatable arrays.
2. The first data-matrix column remains the reference series, matching the R
   package.
3. `QCSIS::qc` is called once with the complete quantile vector instead of
   once per quantile. This avoids repeated sorting while producing the same
   estimates.
4. Simulated correlations are retained only for one wavelet level at a time,
   reducing memory from `O(n_sim * J * n_quantiles)` to
   `O(n_sim * n_quantiles)`.
5. Confidence limits use R's default type-7 sample quantile at probabilities
   0.025 and 0.975.
6. A private RNG object prevents the analysis from depending on a compiler's
   implementation of `random_number`. An optional seed provides deterministic
   tests and repeatable analyses.
7. Plotting, lattice objects, grid drawing, palettes, and R data-frame class
   infrastructure are omitted.

## Dependencies

The attached modern Fortran translations are bundled unchanged as local FPM
path dependencies:

- QCSIS for quantile correlation;
- waveslim for MODWT MRA.
