# Porting notes

## Design

The R package's vectors, lists, S4 distribution objects, and model objects are
represented by explicit Fortran derived types:

- `sn_uv_params`, `st_uv_params`
- `sn_mv_params`, `st_mv_params`
- `sun_params`
- `selm_result`, `grouped_fit_result`

All public numerical routines use `real(dp)`, where `dp = kind(1.0d0)`, and
return explicit status codes when failure is possible.

## Direct translations

The following are formula-level translations from upstream code:

- univariate SN/ESN, ST, and SC densities and distribution functions
- transformation-based random generators
- Owen T, zeta derivatives, cumulants, modes, and centered parameters
- multivariate SN/ST/SC densities and latent-variable random generators
- SUN density algebra and distribution transformations
- Q penalty
- matrix half-vectorization and duplication operators

## Self-contained numerical replacements

The upstream package imports `mnormt`, `numDeriv`, and `quantreg`. This port has
no external numerical dependencies:

- dimensions 1 and 2 use direct normal probability formulas
- higher-dimensional normal probabilities use antithetic Halton integration
- multivariate Student-t probabilities use deterministic Halton integration
- numerical gradients/Hessians use central differences
- likelihood fitting uses Nelder-Mead with transformed scale/degrees-of-freedom
- grouped data are fitted from interval probability differences
- SUN moments use deterministic truncated-normal simulation
- the Galton-Moors initializer uses a deterministic distribution grid
- multivariate fits combine marginal penalized fits with an empirical
  correlation estimate

These replacements preserve the statistical meaning and deterministic behavior,
but not the exact integration/optimization path of the R dependencies.

## Parameter conventions

`omega` and `Omega` are scale parameters, not necessarily the ordinary standard
deviation/covariance after skewing. `alpha` is the direct slant parameter.
For ST fitting, `nu = 1 + exp(theta)` guarantees positive degrees of freedom.

The default SN/ESN/ST/SC regression location is `X beta`. Use
`predict_selm(..., mean_response=...)` when the distributional mean is required
and exists.

## Licensing

Upstream declares `GPL-2 | GPL-3`. Every translated Fortran source therefore
uses the SPDX expression `GPL-2.0-only OR GPL-3.0-only`. The original archive is
included unchanged for provenance.
