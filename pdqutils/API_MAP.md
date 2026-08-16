# API map

| R PDQutils export | Fortran API | Status |
|---|---|---|
| `AS269` | `as269`, `as269_orders`, `as269_vector` | translated |
| `moment2cumulant` | `moment2cumulant` | translated |
| `cumulant2moment` | `cumulant2moment` | translated |
| `dapx_edgeworth` | `dapx_edgeworth`, `dapx_edgeworth_vec` | translated |
| `papx_edgeworth` | `papx_edgeworth`, `papx_edgeworth_vec` | translated |
| `dapx_gca` | `dapx_gca`, `dapx_gca_vec` | translated |
| `papx_gca` | `papx_gca`, `papx_gca_vec` | translated |
| `qapx_cf` | `qapx_cf`, `qapx_cf_vec` | translated |
| `rapx_cf` | `rapx_cf` | translated |

The R `basis` strings are represented by integer constants
`gca_normal`, `gca_gamma`, `gca_beta`, `gca_arcsine`, and `gca_wigner`.
`gca_basis_from_name` is provided for callers that prefer string selection.

R's `basepar` list is represented by optional keyword arguments:

- gamma: `shape=`, `scale=`
- beta: `shape1=`, `shape2=`

R graphics, package documentation machinery, and vignette rendering are not
computational APIs and are not part of the Fortran library.
