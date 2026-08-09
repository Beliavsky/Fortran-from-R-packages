# Translation coverage

## Directly translated numerical behavior

| Upstream feature | Fortran status |
|---|---|
| Ordinary NLS initialization | Implemented with standalone bounded LM |
| Orthogonal residual objective | Implemented |
| Per-observation `optimize()` projection | Implemented with Brent minimization |
| `window` local projection intervals | Implemented |
| `extend` global projection interval | Implemented |
| parameter bounds | Implemented |
| observation weights | Implemented |
| `fixed` parameters | Implemented with upstream ONLS-stage semantics |
| `x0`, `y0` | Implemented |
| vertical/orthogonal residuals and RSS | Implemented |
| numerical model gradient / covariance / SE | Implemented |
| `logLik.onls`, `logLik_o` numerical formulas | Implemented |
| `check_o` orthogonality angles | Implemented without plotting |
| `NIST()` | Chwirut2 retained as an executable example/regression case |

## Deliberate representation changes

R formulas, model frames, environments, parameter-name lookup, self-start
objects, NA actions, and S3 classes are replaced by explicit numeric arrays and
a typed model callback.  This is an API change, not a change in the orthogonal
least-squares objective.

The upstream `minpack.lm::nls.lm` dependency is not bundled in `onls`; this
translation provides a self-contained LM/Gauss-Newton implementation.  Thus
iteration counts and exact convergence messages are not expected to match
MINPACK even when the final solution does.

R's `stats::optimize` uses a Brent-style bounded minimizer.  The Fortran version
uses a standalone Brent implementation with the same default tolerance used by
`onls` (`sqrt(.Machine$double.eps)`).

## Source quirks intentionally preserved

- `fixed` parameters do not constrain the preliminary ordinary NLS fit.  They
  are reset to their supplied starting values before and during ONLS.
- Predictor/response data are sorted before ONLS.  The R source does not reorder
  the precomputed weight vector.  This is the default Fortran behavior; it can
  be disabled with `mimic_r_unsorted_weights=.false.`.
- Stored `resid_o` values are weighted orthogonal distances.  Consequently the
  upstream `logLik_o()` formula applies weights again; the Fortran helper
  preserves that calculation.  Raw distances are additionally exposed as
  `distance_o`.
- The formal `jac` argument of upstream `onls()` is not passed to either
  `nls.lm` call and therefore has no numerical effect; no corresponding unused
  argument is included in the Fortran API.

## Not translated

These are R/model-object features rather than the core fit engine:

- formula parsing, `model.frame`, `getInitial`, self-start model discovery;
- S3 print/plot methods and graphics;
- NA/subset handling and `napredict` reconstruction;
- S3 `predict` wrapper (prediction is simply a model callback evaluation);
- `profile.onls` / profile-based `confint.onls` and their R spline/t/F
  distribution orchestration.

The complete supplied R package is retained under `original/onls-master/` for
provenance and for future extension of profile-based inference if desired.
