# Porting notes

## R objects

The R package returns S3 classes such as `pam`, `agnes`, `diana`, `fanny`, and
`dissimilarity`. The Fortran port uses explicit derived result types instead of
runtime class dispatch, attributes, formulas, or `...` arguments.

## Distances and missing data

Observations are rows. IEEE NaN values are omitted coordinate-by-coordinate in
numeric distances. Mixed data are passed in one real matrix plus an integer
variable-type vector; categorical values are numeric codes.

## Algorithm adaptations

- PAM uses the classical BUILD initialization and exhaustive SWAP search.
- CLARA uses deterministic portable sampling and evaluates candidate medoids on
  the full data set.
- FANNY is implemented as fuzzy c-medoids with the standard membership update.
- AGNES reuses a portable Lance-Williams engine. It supports the five most used
  linkage families; R's flexible and generalized-average parameterizations are
  not included in this release.
- DIANA uses the original splinter-group idea and builds a postorder merge tree.
- MONA implements recursive monothetic binary splitting, but does not reproduce
  every upstream banner/object field.
- `clus_gap` generates references in the original coordinate bounding box. The
  R option that first rotates to scaled principal-component coordinates is not
  included.
- The enclosing ellipsoid uses a Khachiyan iteration. If its conservative
  iteration cap is reached, the result is rescaled to guarantee enclosure and
  returned with a descriptive success message.

## Omitted code

Plotting (`clusplot`, banner plots, tree plots), interactive menus, S3 print and
summary methods, localization catalogs, and R data-object plumbing are omitted.
No computational result depends on plotting code.
