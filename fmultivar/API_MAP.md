# API map

This table maps the original R package routines to the translated Fortran API.

| R routine | Fortran procedure or type | Status |
|---|---|---|
| `dnorm2d` | `dnorm2d` | Direct density formula; tested |
| `pnorm2d` | `pnorm2d` | Conditional adaptive quadrature; tested |
| `rnorm2d` | `rnorm2d` | Direct correlated simulation; tested through RNG moments |
| `dt2d` | `dt2d` | Direct density formula; tested |
| `pt2d` | `pt2d` | Conditional adaptive quadrature; tested |
| `rt2d` | `rt2d` | Scale-mixture simulation; tested through multivariate t RNG paths |
| `dcauchy2d` | `dcauchy2d` | Student-t with one degree of freedom; tested |
| `pcauchy2d` | `pcauchy2d` | Student-t CDF with one degree of freedom; tested |
| `rcauchy2d` | `rcauchy2d` | Student-t RNG with one degree of freedom; tested through compatibility paths |
| `delliptical2d` | `elliptical2d_density` | All seven generator families; tested |
| `.gfunc2d` | Internal branch in `elliptical2d_density` | Direct formulas; tested through all families |
| `dmvnorm` | `dmvnorm`, `mvnorm_pdf`, `mvnorm_logpdf` | Tested |
| `pmvnorm` | `pmvnorm`, `mvnorm_rect_prob` | Reproducible Monte Carlo with standard error; tested |
| `qmvnorm` | `qmvnorm`, `mvnorm_equicoordinate_quantile` | Common-random-number bisection; tested |
| `rmvnorm` | `rmvnorm`, `mvnorm_rng` | Tested |
| `dmvt` | `dmvt`, `mvt_pdf`, `mvt_logpdf` | Tested |
| `pmvt` | `pmvt`, `mvt_rect_prob` | Reproducible Monte Carlo with standard error; tested |
| `qmvt` | `qmvt`, `mvt_equicoordinate_quantile` | Common-random-number bisection; tested |
| `rmvt` | `rmvt`, `mvt_rng` | Tested |
| `dmsn` / `dmvsnorm` | `dmsn`, `mvsnorm_pdf` | Direct skew-Normal formula; tested |
| `pmsn` / `pmvsnorm` | `pmsn`, `mvsnorm_rect_prob` | Reproducible Monte Carlo; tested |
| `rmsn` / `rmvsnorm` | `rmsn`, `mvsnorm_rng` | Latent-selection simulation; tested |
| `dmst` / `dmvst` | `dmst`, `mvst_pdf` | Direct skew-t formula; tested |
| `pmst` / `pmvst` | `pmst`, `mvst_rect_prob` | Reproducible Monte Carlo; tested |
| `rmst` / `rmvst` | `rmst`, `mvst_rng` | Latent scale-mixture simulation; tested |
| `dmsc` | `dmsc` | Skew-t with `nu=1`; tested |
| `pmsc` | `pmsc` | Skew-t probability with `nu=1`; tested |
| `rmsc` | `rmsc` | Skew-t simulation with `nu=1`; tested |
| `.mnFit` | `fit_multivariate_normal` | Sample mean/covariance fit; tested |
| `msnFit` | `fit_skew_normal`, `msn_fit` | Numerical MLE analogue; tested in 1-D and 2-D |
| `mstFit` | `fit_skew_t`, `mst_fit` | Fixed and free `nu`; tested, including 2-D fixed-`nu` fit |
| `mscFit` | `fit_skew_cauchy`, `msc_fit` | Numerical MLE analogue; tested in 1-D and 2-D |
| `mvFit` | `mv_fit`, `mvfit` | Dispatch wrapper; tested |
| `adapt` | `adapt_integrate_nd`, `adapt` | Halton quasi-Monte Carlo analogue, dimensions 1-20; tested in 2-D and 3-D |
| `integrate2d` | `integrate2d_rule`, `integrate2d` | Original unit-square nine-point rule; tested |
| none | `adapt_integrate2d` | Additional adaptive 2-D Gauss-Legendre routine; tested |
| `grid2d` | `grid2d`, type `grid_coordinates` | Tested |
| `gridData` | `make_grid_data`, `griddata`, type `grid_data` | Numeric container only; tested |
| `density2d` | `density2d` | Axis-aligned Normal KDE; tested by numerical integration |
| `hist2d` | `hist2d` | Tested by count conservation |
| `squareBinning` | `square_binning`, `squarebinning`, type `binning_result` | Tested by count conservation |
| `hexBinning` | `hex_binning`, `hexbinning`, type `binning_result` | Tested by count conservation |
| plot and slider methods | none | Excluded as plotting/R infrastructure |

The Fortran fit result replaces the R `fDISTFIT` object with the plain derived
type `skew_fit_result`.
