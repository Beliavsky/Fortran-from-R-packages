# API coverage

Coverage basis: exported computational R functions.

Package status: **substantial**.

Mapped functions: **8 of 8 (100.0%)**.

The percentage records whether each exported computational operation has a
meaningful public Fortran implementation. It does not mean that R S4 dispatch,
object layouts, warnings, plotting, callbacks, or every optional argument are
identical.

| R function | Fortran procedure(s) | Status | Main differences |
| --- | --- | --- | --- |
| `kepdf` | `kepdf_fit`, `kepdf_estimate` | substantial | Numeric KDE translated; S4/data-frame interface omitted. |
| `h.norm` | `h_norm` | complete | Direct numerical translation. |
| `hprop2f` | `hprop2f` | complete | Direct pilot-density/adaptive-bandwidth translation. |
| `pdfCluster` | `pdfcluster_fit`, `pdfcluster_regraph_pairs`, `density_profile`, `core_labels_from_profile` | substantial | Numeric graph/profile/core algorithm translated; S4/dendrogram/Qhull-option behavior omitted. |
| `pdfClassification` | `pdf_classification` | substantial | Staged classifier translated; zero replaces R `NA`; typed API replaces R object mutation. |
| `groups` | `groups_at_stage` | substantial | Numeric stage extraction translated; zero replaces R `NA`. |
| `dbs` | `dbs_from_labels`, `dbs_from_pdfcluster` | substantial | Density-based silhouette translated; arbitrary R bandwidth-function lookup and S4 output omitted. |
| `adj.rand.index` | `adjusted_rand_index` | complete | Direct contingency-table formula; degenerate zero denominator returns IEEE NaN. |

Excluded from the denominator are plotting methods, `summary`/`show` methods,
class declarations, and presentation-only interfaces.
