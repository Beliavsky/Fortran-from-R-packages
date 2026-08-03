# Porting notes

## R objects

The R package returns data frames or lists and uses S3 clustering objects. The
Fortran translation returns a typed `portfolio_result`; asset names are left to
the calling application.

## Clustering dependencies

The upstream package imports `fastcluster` and `cluster`. This project is
self-contained and adapts the portable Lance-Williams, DIANA, cut-tree, and
gap-statistic components developed in the earlier Fortran translations of
those packages. Supported agglomerative linkages are the four methods used by
HierPortfolios: single, complete, average, and Ward.D2.

The common clustering engine uses `O(n^2)` memory and `O(n^3)` worst-case time.
This is slower than the specialized upstream `fastcluster` kernels for large
asset universes but is simple, deterministic, and portable.

## Gap statistic

The R package delegates to `cluster::clusGap`, whose default reference space
uses scaled PCA. This translation uses deterministic uniform simulations over
the axis-aligned bounding box of the correlation-distance features. It applies
the Tibshirani one-standard-error selection rule and then enforces at least two
clusters, matching the portfolio functions' intent. Set `clusters` explicitly
when exact cluster-count control is required.

## HERC

The upstream implementation stores risk values across every cut of the tree.
The Fortran implementation computes the equivalent allocation more directly:
terminal clusters receive inverse-variance portfolios, and each successive
hierarchical split allocates capital inversely to the aggregate risk of its two
child groups.

## DHRP bounds

The upstream R code applies aggregate lower and upper bounds during each split.
That iteration can be numerically fragile and does not always guarantee final
per-asset feasibility. The Fortran routine preserves the DIANA and `tau` split
logic, then projects the unconstrained weights onto the box-constrained simplex
`lb <= w <= ub`, `sum(w)=1`. This guarantees the documented asset bounds when
they are feasible.

## Plotting and data

Dendrogram plotting is omitted. The original `.RData` files are retained under
`original/` but are not parsed by the compiled library.
