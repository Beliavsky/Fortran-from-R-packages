# API coverage

This document describes computational coverage relative to upstream `spdep`
1.4-2. The Fortran API intentionally uses arrays and derived types rather than
R S3 classes, formula objects, data frames, `sf` objects, or R method dispatch.

## Implemented computational areas

### Neighbor graphs and geometry

- `cell2nb` rook/queen grids, with optional toroidal wrapping
- `dnearneigh` planar and upstream-compatible WGS84 distance bands
- `knearneigh`, `knn2nb`
- `nbdists`
- Gabriel graph and relative-neighborhood graph
- Delaunay-neighbor graph (`tri2nb`) using an empty-circumcircle algorithm
- cardinalities (`card`), weak connected components, reciprocal-symmetry tests
- symmetrization, self-link inclusion/removal
- union, intersection, set difference, and upstream-style graph complement
- directed link addition/removal
- exact graph-lag neighbor lists, cumulative lag unions, and all-pairs unweighted graph distances
- block-level neighbor expansion (`nb2blocknb`) and 2D coordinate rotation (`Rotation`)
- dense logical adjacency conversion

### Spatial weights

- `nb2listw` styles `B`, `W`, `C`, `U`, `S`, and `minmax`
- optional general link weights (`glist` analogue)
- dense `nb2mat`, `mat2listw`, `listw2mat`
- vector and matrix spatial lags
- `S0`, `S1`, `S2` spatial-weight constants and `Szero`
- symmetric `listw2U`
- distance-decay `nb2listwdist` weights (`idw`, `exp`, `dpd`)
- `autocov_dist` distance-weighted spatial autocovariates

### Spatial statistics

- global Moran's I
- analytical Moran test under randomization or normality moments
- deterministic seeded Moran Monte Carlo test
- bivariate Moran statistic
- local Moran statistic with conditional or unconditional moments
- local bivariate Moran statistic
- global Geary's C
- Local Geary `local_geary` (upstream `localC`) deterministic contributions
- analytical Geary test under randomization or normality moments
- deterministic seeded Geary Monte Carlo test
- global Getis-Ord G analytical test
- local Getis-Ord G and G-star, selected by self-inclusion in the weights
- global and local Lee spatial association
- LOSH, including the upstream chi-square inference when `a=2`
- univariate binary join-count analytical test under fixed or free sampling
- global and local empirical-Bayes rate smoothing
- Choynowski Poisson tail probabilities
- empirical-Bayes-adjusted Moran statistic

### Regionalization

- neighbor feature costs (`nbcosts`) using Euclidean, Manhattan, maximum,
  Canberra, binary, or Minkowski distances
- within-group dispersion (`ssw`)
- Prim minimum spanning tree (`mstree`)
- edge pruning gains (`prunecost`)
- iterative SKATER-style tree pruning (`skater_groups`) with an optional
  minimum group size

### Weighted multivariate spatial autocorrelation

- `linearised_diffusive_weights`
- `metropolis_hastings_weights`
- `iterative_proportional_fitting_weights`
- `graph_distance_weights`
- global `spatialdelta` normal approximation including skewness and excess
  kurtosis moments
- `localdelta`
- Cornish-Fisher correction

## Deliberately not translated

The following are outside this package's Fortran computational API or require
R-specific object systems/external spatial geometry infrastructure:

- plotting, printing, summaries whose only role is presentation, palettes, and
  interactive helpers
- `sf`, `s2`, `sp`, and polygon-object import/conversion methods, including
  polygon contiguity construction from R geometry objects
- R formula/model-frame methods, `lm` diagnostics, and tests tied to R model
  objects; upstream notes that most spatial regression fitting moved to
  `spatialreg`
- serialization, row names, factor handling, NA-action dispatch, R options,
  namespace plumbing, and S3 replacement/subsetting methods
- exact and saddlepoint tests that are tightly coupled to R linear-model and
  eigen-object workflows
- optional plotting-coordinate/scree routines from the spatial-delta API
- local permutation-reporting layers such as `localC_perm`, `localG_perm`, and
  `localmoran_perm`; the deterministic local statistics are translated
- multicolour/local categorical reporting layers that primarily build R tables
  around repeated permutation kernels

These omissions are explicit rather than silent. The translated API focuses on
portable numerical kernels that are useful directly from Fortran and can later
be wrapped from Python, R, or other languages.
