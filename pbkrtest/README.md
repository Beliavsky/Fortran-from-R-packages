# pbkrtest

Modern free-form Fortran translation of the portable computational kernels in
R package **pbkrtest 0.5.5**.

The translated API focuses on the numerical work behind Kenward-Roger and
Satterthwaite inference, parametric-bootstrap calibration, likelihood-ratio
calibration, covariance-component decomposition, and restriction/model-matrix
linear algebra. R formulas, S3 methods, model refitting, data-frame output,
printing, plotting, parallel R clusters, and R-specific RNG wrappers are out of
scope.

## Repository layout and dependencies

Place this directory at the root of `Fortran-from-R-packages`, alongside:

```text
pbkrtest/
rfortran-core/
rfortran-linalg/
numDeriv/
```

`fpm.toml` uses sibling path dependencies and does not vendor their source.
`rfortran-linalg` supplies the shared pinned `fortran-lapack` dependency, so the
pbkrtest package itself contains no `-lblas`, `-llapack`, or system-library
link directives.

## Build and test

With FPM and gfortran available:

```text
fpm build
fpm test
fpm run --example basic_example
```

The package is free-form Fortran and uses a single real kind, `dp`, re-exported
through the public `pbkrtest` module from `rfortran-core`'s `r_kinds` module.

## Main public API

- `build_sigma_g`: construct `Sigma`, covariance derivative matrices `G`, and
  covariance parameters `gamma` from random-effect design blocks.
- `vcov_adjust_kr`: Kenward-Roger adjusted fixed-effect covariance and its
  `P`, `W`, and information matrices.
- `kr_adjust`: adjusted F statistic, scaling, and denominator degrees of
  freedom for a restriction matrix.
- `lb_ddf`, `ddf_lb_scalar`: Kenward-Roger denominator-df helpers.
- `compute_auxiliary_numeric`: callback-based numerical Hessian and covariance
  Jacobian calculation using the shared `numDeriv` translation.
- `satterthwaite_test`, `get_fstat_ddf`: Satterthwaite F approximation.
- `likelihood_ratio_test`: LRT statistic and chi-square tail probability.
- `bootstrap_p_values`: direct bootstrap p-value plus Bartlett, gamma, and
  moment-matched F calibrations from an already-generated reference sample.
- `compare_column_space`, `orthogonal_complement`, `force_full_rank`,
  `make_model_matrix`, `make_restriction_matrix`: linear-algebra utilities used
  to translate model nesting into restrictions and back.

See `API_COVERAGE.md` for a mapping from upstream R routines to Fortran APIs and
for explicitly omitted interfaces.

## Design notes

`compute_auxiliary_numeric` deliberately receives deviance and vectorized
fixed-effect covariance callbacks. This preserves the numerical algorithm in
`compute_auxiliary()` without embedding R's `lmerMod` environment machinery.
Similarly, `bootstrap_p_values` translates the calibration step in
`PBcompute_p_values()` but not R/lme4 model simulation or parallel refitting.

Every maintained procedure dummy argument has explicit `INTENT` or `VALUE` and
a trailing FORD documentation comment. NaN construction/tests use the intrinsic
IEEE arithmetic module; no fast-math assumptions are required or enabled.

## License and provenance

The upstream package declares `GPL (>= 2)`. This translation is therefore GPL-2.0-or-later.
See `LICENSE`, `COPYING.GPL-2`, `COPYING.GPL-3`, and `NOTICE.md`.
