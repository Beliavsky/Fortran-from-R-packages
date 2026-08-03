# Translation coverage

## Direct computational exports represented

| R export | Fortran representation |
|---|---|
| `agnes` | `agnes`, `agnes_distance` |
| `clara` | `clara` |
| `daisy` | `daisy`, `daisy_mixed` |
| `diana` | `diana`, `diana_distance` |
| `fanny` | `fanny` |
| `mona` | `mona` |
| `pam` | `pam`, `pam_distance` |
| `silhouette` | `silhouette` |
| `sortSilhouette` | `sort_silhouette` |
| `clusGap` | `clus_gap` |
| `maxSE` | `max_se` |
| `ellipsoidhull` | `ellipsoidhull` |
| `ellipsoidPoints` | `ellipsoid_points` |
| `predict.ellipsoid` | `predict_ellipsoid` |
| `volume.ellipsoid` | `volume_ellipsoid` |
| `medoids` | `medoids` |
| `coef.hclust`, `coefHier` | stored coefficient and `coef_hier` |
| `meanabsdev` | `meanabsdev` |
| `sizeDiss` | `size_diss` |
| triangular-index helpers | `lower_to_upper_tri_inds`, `upper_to_lower_tri_inds` |

## Omitted exports

`bannerplot`, `clusplot`, and `pltree` are graphical. R print, plot, summary,
menu, and object-conversion methods are also omitted.

## Partial compatibility

AGNES flexible/generalized-average linkage and scaled-PCA gap-reference space are
not included. See `PORTING.md` for numerical differences.
