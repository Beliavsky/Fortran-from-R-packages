# Omitted R-specific functionality

The original package contains no standalone plotting implementation in its `R/`
or `src/` directories. Plotting calls appearing in examples and vignettes depend
on external R packages such as `igraph`, `corrplot`, and `viridis`; those calls
are not part of this Fortran translation.

The following infrastructure is also intentionally omitted:

- Rcpp registration and generated wrappers;
- roxygen documentation objects and R package namespace machinery;
- progress bars and S3/list presentation conventions;
- benchmark, website, vignette rendering, and CI configuration.

All exported numerical graph operators, metrics, and graph-learning routines
from `NAMESPACE` are represented in the Fortran API. The non-exported
`learn_cospectral_graph` computational routine is included as well.
