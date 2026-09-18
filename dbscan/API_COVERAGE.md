# API coverage

Coverage basis: exported computational R functions from the supplied `dbscan`
1.2.6 `NAMESPACE`.  Presentation-only exports (`augment`, `clplot`, `glance`,
`hullplot`, `kNNdistplot`, and `tidy`) are excluded.  Internal helpers and
unexported S3 presentation methods are also excluded.

**Status: substantial — 25 of 25 (100%) functions mapped.**

The 100% mapping fraction is not a claim of full R compatibility.  The largest
remaining computational differences are:

- nearest-neighbor searches are exact brute force rather than ANN kd-tree or
  priority/approximate search;
- `hdbscan` and `glosh` implement the upstream condensed-tree stability,
  GLOSH, and `cluster_selection_epsilon` computations for numeric matrix input,
  but R `dist` input and optional presentation-tree objects are omitted;
- `extractFOSC` supports unsupervised extraction, symmetric full-matrix
  pairwise constraints, mixed `alpha` objectives, and unstable-branch pruning;
  R list/dist-vector constraint forms, automatic constraint repair, and hclust
  object augmentation are omitted;
- `as.dendrogram`/`as.reachability` provide numeric translated structures, not
  arbitrary R dendrogram/hclust/S3 object semantics;
- R `dist` class metadata, recycling, partial argument matching, warnings, and
  S3 dispatch are generally omitted.

The authoritative per-function mapping, statuses, source locations, and notes
are stored in `fpm.toml` and mirrored in the README table.
