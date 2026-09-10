# API coverage

Coverage basis: distinct exported computational R functions in upstream
`NAMESPACE`. Printing and package-documentation methods are excluded.

| R function | Fortran mapping | Status | Main differences |
| --- | --- | --- | --- |
| `biv_kld` | `bivkld_kernel:biv_kld` | substantial | Numeric n-by-2 arrays replace R matrices/data frames; R exceptions become `info` plus NaN; SCV uses a streamlined pilot/optimizer; S3 construction is replaced by optional `bivkld_estimate`. |
| `biv_kld_matrix` | `bivkld_kernel:biv_kld_matrix` | substantial | `biv_sample(:)` replaces R list/group splitting; optional H supplies all group bandwidths together; names/attributes/S3 class are omitted; numerical common-grid algorithm is retained. |
| `biv_kld_discrete` | `bivkld_exact:biv_kld_discrete_vector`, `biv_kld_discrete_matrix` | substantial | Vectors and rank-2 tables are supported; higher-rank R arrays are not a separate overload; invalid inputs return NaN instead of raising an R error. |
| `biv_kld_normal` | `bivkld_exact:biv_kld_normal` | complete | Explicit two-by-two algebra replaces `solve`/`determinant`; invalid inputs return NaN. |
| `biv_kld_pareto2` | `bivkld_exact:biv_kld_pareto2` | complete | Formula preserved; invalid inputs return NaN. |
| `biv_kld_independent_weibull` | `bivkld_exact:biv_kld_independent_weibull` | complete | Formula preserved; invalid inputs return NaN. |

The coverage count is 6/6 = 1.0. It measures whether each exported
computational R function has a meaningful public Fortran implementation, not
whether all R object-system semantics and dependency options are reproduced.
