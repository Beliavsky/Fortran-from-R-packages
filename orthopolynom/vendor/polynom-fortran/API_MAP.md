# API map

| Upstream R API | Fortran API | Notes |
|---|---|---|
| `polynomial()` | `polynomial(...)` | Generic constructor for real/integer scalars and vectors; no-argument form returns `x`. |
| `as.polynomial()` / `is.polynomial()` | `type(polynomial_t)` | Static typing removes runtime coercion/class checks. |
| `polylist()` / `as.polylist()` / `is.polylist()` | `type(polylist_t)` | Uses an allocatable array of `polynomial_t`. |
| `Ops.polynomial` | overloaded `+ - * / ** == /=`; `poly_rem`, `poly_divmod` | Fortran has no `%%` operator, so remainder is named. |
| `Summary.polynomial`, `Summary.polylist` | `sum_polynomials`, `product_polynomials` | Explicit reductions. |
| `Math.polynomial` | `round_coefficients`, `significant_coefficients`, `floor_coefficients`, `ceiling_coefficients`, `truncate_coefficients` | Explicit numerical transformations. |
| `coef.polynomial` | `%coefficients()`, `coef_polynomial` | Returns a copy of the coefficient vector. |
| `predict.polynomial` | `%evaluate()`, `predict_polynomial` | Scalar and vector overloads. |
| `as.function.polynomial` | `%evaluate()` | Procedure closure generation is unnecessary in Fortran. |
| `as.character.polynomial` | `%to_string()` | Formatting only; no R print method. |
| `deriv.polynomial` | `derivative`, `deriv_polynomial` | Optional derivative order supported. |
| `deriv.polylist` | `derivative_polylist` | Elementwise mapping. |
| `integral.polynomial` | `integral_polynomial`, `definite_integral`, generic `integral` | Antiderivative or two-limit integral. |
| `integral.polylist` | `integral_polylist` | Elementwise mapping. |
| `change.origin` | `change_origin` | Binomial-transform implementation. |
| `monic` | `monic` | Returns a status for the zero polynomial. |
| `poly.calc(x)` | `poly_calc(x)`, `poly_from_roots` | Real-root construction. |
| `poly.calc(x,y)` | `poly_calc(x,y)`, `poly_from_values` | Lagrange interpolation. |
| matrix `y` in `poly.calc` | `poly_from_values_matrix` | One polynomial per matrix column. |
| `poly.from.zeros`, `poly.from.roots` | `poly_from_zeros`, `poly_from_roots` | Aliases retained. |
| `poly.from.values` | `poly_from_values` | Direct mapping. |
| `poly.orth` | `poly_orth`, `orthogonal_polynomials` | Same three-term recurrence and normalization. |
| `solve.polynomial` | `solve_polynomial`, `polynomial_roots` | Aberth iteration replaces companion-matrix eigenvalues. |
| `summary.polynomial` | `summarize_polynomial` | Typed roots/stationary/inflexion result. |
| `GCD`, `LCM` | generic `gcd`, `lcm`; `polynomial_gcd`, `polynomial_lcm`, list variants | Pair and polylist forms. |
| `[`, `c`, `rep`, `unique` for polylist | native Fortran array operations | No wrapper needed. |
| plotting methods | omitted | Non-computational graphics. |
