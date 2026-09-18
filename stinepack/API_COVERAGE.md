# API coverage

Coverage basis: distinct exported computational R functions in `stinepack` 1.5.
Presentation-only functions and internal helpers are excluded. The R function
`na.stinterp` is a dispatch-only S3 generic and is therefore excluded; the
registered computational method `na.stinterp.default` is included.

Package status: **substantial**  
Mapped functions: **4 of 4 (100%)**

| R function | R source | Fortran mapping | Status | Notes |
|---|---|---|---|---|
| `na.stinterp.default` | `upstream/R/na.stinterp.R` | `stinepack_api: na_stinterp` | substantial | Rank-1/rank-2 numerical gap filling, default/equally spaced or explicit coordinates, `na_rm`, and method selection are implemented. S3/class/time attributes and arbitrary `...` forwarding are omitted. |
| `parabolaSlopes` | `upstream/R/parabolaSlopes.R` | `stinepack_api: parabola_slopes` | complete | Direct numerical formula translation for valid knot vectors. |
| `stinemanSlopes` | `upstream/R/stinemanSlopes.R` | `stinepack_api: stineman_slopes` | complete | Both scaled and unscaled formulas are translated. |
| `stinterp` | `upstream/R/stinterp.R` | `stinepack_api: stinterp` | substantial | Core rational interpolation, all three slope choices, explicit slopes, endpoint tolerance, and out-of-range NaNs are implemented. R errors/list semantics are represented with a derived result type and status/message. |

The 100% fraction means that every function in the stated computational
coverage basis has a meaningful mapping. It does not imply full R-interface or
R-object compatibility.
