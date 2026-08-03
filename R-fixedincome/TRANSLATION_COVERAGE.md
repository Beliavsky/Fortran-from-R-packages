# Translation coverage

The upstream `NAMESPACE` contains 41 directly exported names. The single
plotting export, `ggspotratecurveplot`, is intentionally omitted. The other 40
exports are represented as follows.

| Upstream export | Fortran representation |
|---|---|
| `interpolation<-` | `set_interpolation` |
| compounding classes | integer constants and `compounding` |
| `as.forwardrate` | `as_forwardrate`, `forwardrate_from_curve` |
| `as.spotrate` | `as_spotrate` |
| `as.spotratecurve` | `as_spotratecurve` |
| `as.term` | `parse_term` |
| `closest`, `first`, `last` | same names |
| `compound`, `discount`, `implied_rate` | same scalar/vector kernels plus typed rate wrappers |
| `compounding` | same name |
| `daycount`, `dib` | same names |
| `fit_interpolation` | same name |
| `forwardrate` | constructor plus curve conversion routines |
| six nonparametric interpolation constructors | same names with dots removed where applicable |
| Nelson-Siegel constructors | same names with dots removed where applicable |
| `interpolate`, `interpolation`, `interpolation_error` | same names |
| `maturities` | same name; returns ordinal dates |
| `parameters` | same name |
| `prepare_interpolation` | validation-oriented typed equivalent |
| `shift` | same name |
| `spotrate`, `spotratecurve`, `term` | same names |
| `todays`, `tomonths`, `toyears` | same names |
| `ggspotratecurveplot` | omitted plotting code |

Non-exported numerical functions for Nelson-Siegel/Svensson evaluation and
objectives were also translated. R-only methods for indexing, formatting,
plotting, and data frames were replaced by ordinary Fortran array access or
explicit procedures rather than duplicated as artificial wrappers.
