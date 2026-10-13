# API coverage and compatibility

Status: **substantial**. Coverage basis: exported computational R functions.
Mapped **15 of 15 (100%)**; untranslated list is empty. This counts mapped
functions, including exported aliases, rather than complete R compatibility.
The manifest and README contain the same 15 entries.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `bSpline` | `splines2` | `b_spline` | substantial |
| `bsp` | `splines2` | `bsp` | substantial |
| `dbs` | `splines2` | `dbs` | substantial |
| `ibs` | `splines2` | `ibs` | substantial |
| `mSpline` | `splines2` | `m_spline` | substantial |
| `msp` | `splines2` | `msp` | substantial |
| `iSpline` | `splines2` | `i_spline` | substantial |
| `isp` | `splines2` | `isp` | substantial |
| `cSpline` | `splines2` | `c_spline` | substantial |
| `csp` | `splines2` | `csp` | substantial |
| `naturalSpline` | `splines2` | `natural_spline` | substantial |
| `nsp` | `splines2` | `nsp` | substantial |
| `nsk` | `splines2` | `nsk` | substantial |
| `bernsteinPoly` | `splines2` | `bernstein_poly` | substantial |
| `bpoly` | `splines2` | `bpoly` | substantial |

## Inventory rules

`upstream/NAMESPACE` exports exactly the above 15 names. The installed pinned
package's `getNamespaceExports("splines2")` returns the same list. Although R
sources use roxygen `@export` for S3 registration, their generated NAMESPACE
has only `S3method(...)` entries for those methods; they are not exported
symbols. They do not enter the explicitly exported computational denominator.
Internal Rcpp engine wrappers and helpers are not exported. Plot/print methods
are presentation-only. Formula dispatch, extraction, and makepredictcall object
plumbing do not need a Fortran emulation.

## Material differences

1. Constructors return a typed matrix/settings/status object rather than an R
   matrix with attributes. No R warnings, class system, attributes, dimension
   names, coercion, recycling, formula parsing, or model-frame processing.
2. NaN rows are preserved; infinite coordinates/knots are rejected explicitly.
   Default boundaries cannot be inferred from empty/all-missing/constant x.
3. Natural trim is restricted to [0,0.5); when boundaries are explicitly supplied,
   trim is ignored. Excessive interior multiplicity above degree+1 is rejected.
4. The upstream 0.5.4 call `cSpline(..., intercept=FALSE, derivs>0, scale=TRUE)`
   fails with `each_row(): incompatible size` because a full-width scale vector
   is applied to a reduced derivative matrix. This translation deliberately
   computes/scales the complete derivative, then drops its first column. The
   corresponding R fixtures use that mathematically equivalent full-basis path.
5. `derivative_spline` adds derivative order to the stored object and allows zero
   as an identity operation. `predict_spline` reconstructs from stored settings.
   These useful APIs do not claim additional exported-R coverage.
6. C-splines have no integral option upstream and none is claimed. Their
   unscaled basis is itself the integral of I-splines. Periodic options exist
   only for B- and M-splines, as in the upstream exported constructors.
7. Upstream's zero-left natural/periodic integral convention is retained,
   even though it differs from integrating a continued basis to a negative x.
8. Degree-zero/one periodic bases with no interior knots are mathematically
   supported. Upstream's df-based construction can throw an Armadillo `subvec`
   error for this case; the translation returns the constant periodic basis
   (or rejects removal of its sole column). The upstream explicit periodic
   M-spline guard `df >= degree` is retained.
9. Very high degree, tightly clustered knots, extreme extrapolation, and poorly
   conditioned cardinal transforms are subject to floating-point conditioning.
   The sibling solver returns a singularity status rather than silently
   producing a cardinal transform after a rejected pivot.

## Implementation map

`src/splines2.f90` exposes every mapped procedure and settings object.
`src/splines2_math.f90` implements spline-specific normalization, cumulative
integrals, closed-form natural transforms, and periodic knot extension.
Cox-de Boor evaluation/derivatives and the cardinal linear solve are reused
through sibling `splines`; sorting/type-7 quantiles through `rfortran-core`.
The shared `splines` natural basis has a different transform and is not used
as a substitute for the splines2 natural basis.
