# API map

## Covariance and location

| R package routine | Fortran routine | Notes |
|---|---|---|
| `Cov`, `CovClassic` | `cov_classic` | Classical mean, covariance, distances |
| `CovMcd` | `cov_mcd` | Random starts, concentration steps, consistency scaling, reweighting |
| `CovMve` / `r_fast_mve` | `cov_mve` | Random elemental starts and minimum-volume subset search |
| `CovOgk` / `covOPW` | `cov_ogk` | Iterated orthogonalized GK; MAD/tau and GK/quadrant options |
| `CovMest`, `covMest` | `cov_mest` | Iteratively reweighted multivariate M estimate |
| `CovSest` / native `sest` | `cov_sest` | High-breakdown biweight S-style estimate |
| `CovMMest` | `cov_mmest` | S initialization followed by efficient redescending M step |
| `CovSde`, `donostah`, native `rlds` | `cov_sde` | Random-direction Stahel-Donoho outlyingness and reweighting |
| `CovMrcd`, `detmrcd` | `cov_mrcd` | Regularized concentration steps with robust diagonal target |
| `covMWcd`, native `fsada` role | `cov_mwcd` | Robust grouped centers and common residual covariance |
| `adjOutlyingness` | `adjusted_outlyingness` | Projection outlyingness with medcouple skew adjustment |
| `mcVT` | `medcouple` | Direct O(n^2) kernel implementation |

All covariance routines return `type(covariance_result)`.

## PCA

| R package routine | Fortran routine |
|---|---|
| `PcaClassic` | `pca_classic` |
| `PcaCov` | `pca_cov` |
| `PcaLocantore` | `pca_locantore` |
| `PcaGrid` | `pca_grid` |
| `PcaProj` | `pca_proj` |
| `PcaHubert` | `pca_hubert` |
| `pca.distances` | `pca_distances` |

PCA routines return `type(pca_result)` with center, scale, loadings,
eigenvalues, scores, score distances, and orthogonal distances.

## Discriminant analysis

| R package routine | Fortran routine |
|---|---|
| `LdaClassic` | `lda_classic_fit` |
| `Linda` | `linda_fit` |
| `LdaPP` | `lda_pp_fit` |
| robust covariance LDA path | `lda_cov_fit` |
| `predict(Lda)` | `lda_predict` |
| `QdaClassic` | `qda_classic_fit` |
| `QdaCov` | `qda_cov_fit` |
| `predict(Qda)` | `qda_predict` |
| `mtxconfusion` | `confusion_matrix` |

Models are stored in `type(lda_model)` and `type(qda_model)`.

## Tests and utilities

| R package routine | Fortran routine |
|---|---|
| `T2.test` | `hotelling_t2_one_sample`, `hotelling_t2_two_sample` |
| `Wilks.test` | `wilks_test` |
| `sqrtm` | `matrix_sqrt` |
| `vecnorm` | `vecnorm` |
| `transform_ilr` | `ilr_transform` |
| covariance flags | `outlier_flags` |
| `cov2cor` use | `covariance_to_correlation` |

Probability values are computed by self-contained regularized gamma and beta
functions.
