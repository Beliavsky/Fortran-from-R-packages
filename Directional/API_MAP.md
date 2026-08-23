# API map

| R function | Fortran API | Status |
|---|---|---|
| `euclid` | `euclid` | translated |
| `euclid.inv` | `euclid_inv` | translated |
| `rotation` | `rotation_matrix` | translated |
| `haversine.dist` | `haversine_dist` | translated |
| `dvm` | `dvm` | translated |
| `dmmvm` | `dmmvm` | translated |
| `dvmf` | `dvmf` | translated |
| `dcardio` | `dcardio` | translated |
| `dwrapcauchy` | `dwrapcauchy` | translated |
| `dwrapnormal` | `dwrapnormal` | translated |
| `dcircbeta` | `dcircbeta` | translated |
| `dcircexp` | `dcircexp` | translated |
| `dcircpurka` | `dcircpurka` | translated |
| `dspcauchy` | `dspcauchy` | translated |
| `dpkbd` | `dpkbd` | translated |
| `dpurka` | `dpurka` | translated |
| `iagd` | `iagd` | translated |
| `desag` (3-D) | `desag3` | translated |
| `pwrapcauchy` | `pwrapcauchy` | translated |
| `pvm` | `pvm_numeric` | numerical integration |
| `rvmf` | `rvmf` | translated |
| `rvonmises` | `rvonmises` | translated |
| `rspcauchy` | `rspcauchy` | translated |
| `circ.summary` | `circ_summary` | core fields translated |
| `circ.cor2` | `circ_cor2` | coefficient translated; p-value placeholder |
| `spher.cor` | `spher_cor_rsq` | diagnostic approximation |
| `mediandir` | `median_direction` | medoid implementation |


## Added in v0.2.0

| R function | Fortran API | Status |
|---|---|---|
| `cardio.mle` | `cardio_mle` | translated; bounded numerical maximization |
| `circexp.mle` | `circexp_mle` | translated; golden-section maximization |
| `vmf.mle` / Rfast equivalent | `vmf_mle` | translated |
| `spcauchy.mle` | `spcauchy_mle` | translated Newton iteration |
| `pkbd.mle` | `pkbd_mle` | translated likelihood with stable numerical gradient search |
| `mixvmf.mle` | `mixvmf_mle` | EM translation with deterministic directional initialization |
| `dirknn` | `dirknn` | translated brute-force directional k-NN |
| `dirda(..., type="vmf")` | `dirda_vmf` | translated |
| `rayleigh` | `rayleigh_test` | asymptotic and Monte Carlo variants |
| `kuiper` | `kuiper_test` | asymptotic and Monte Carlo variants |
| `fb.saddle` | `fb_saddle` | translated saddlepoint approximation |
| `dkent` | `dkent` | translated |
| `rbingham` | `rbingham` | translated ACG rejection sampler |
| `rkent` | `rkent` | compatible rejection sampler; simpler proposal than upstream |
| `matrixfisher.mle` | `matrixfisher_mle` | translated sample-mean SVD |
| `rmatrixfisher` | `rmatrixfisher` | translated via Bingham/quaternion construction |


## Added in v0.3.0

| R function | Fortran API | Status |
|---|---|---|
| `dESAGd` / `desag` | `desag` | arbitrary-dimensional ESAG density |
| `rESAGd` / `resag` | `resag` | arbitrary-dimensional ESAG simulation |
| `.parameter` | `esag_parameters` | translated ESAG gamma-to-eigenvalue/rotation map |
| `rpkbd` | `rpkbd` | translated rejection sampler |
| `dmixspcauchy` | `dmixspcauchy` | translated |
| `dmixpkbd` | `dmixpkbd` | translated |
| `rmixspcauchy` | `rmixspcauchy` | translated; preserves component IDs |
| `rmixpkbd` | `rmixpkbd` | translated; preserves component IDs |
| `mixspcauchy.mle` | `mixspcauchy_mle` | EM with weighted component updates |
| `mixpkbd.mle` | `mixpkbd_mle` | EM with weighted PKBD component updates |
| `kent.mle` | `kent_mle` | Kent orientation plus constrained numerical kappa/beta fit |
| `circ.dcor` | `circ_dcor` | translated through generic distance correlation |
| `spher.dcor` | `spher_dcor` | translated through generic distance correlation |
| `Rfast::dcor` dependency role | `distance_cor` | standalone biased sample distance correlation |
| `knn.reg` | `knn_reg` | scalar/vector response; Euclidean predictors; arithmetic/harmonic response aggregation |
| `knnreg.tune` | `knnreg_tune` | deterministic K-fold tuning; Euclidean or spherical response loss |


## Added in v0.4.0

| R function | Fortran API | Status |
|---|---|---|
| `iag.mle` | `iag_mle` | translated for general dimension; direct IAG likelihood optimization |
| `sipc.mle` | `sipc_mle` | translated likelihood; deterministic pattern-search optimizer |
| `cipc.mle` | `cipc_mle` | translated Newton iteration |
| `gcpc.mle` | `gcpc_mle` | translated GCPC likelihood with positive-rho parameterization |
| `embed.perm` | `embed_perm` | translated two-sample permutation test |
| `hcf.perm` | `hcf_perm` | translated two-sample high-concentration permutation test |
| `het.perm` | `het_perm` | translated heterogeneous vMF permutation test |

## Added in v0.5.0

| R function | Fortran API | Status |
|---|---|---|
| `iag.reg` | `iag_reg` | translated; self-contained likelihood optimization |
| `sipc.reg` | `sipc_reg` | translated; self-contained likelihood optimization |
| `cipc.reg` | `cipc_reg` | translated through the circular spherical-Cauchy regression likelihood |
| `spcauchy.reg` | `spcauchy_reg` | translated for arbitrary directional dimension |
| `pkbd.reg` | `pkbd_reg` | translated for arbitrary directional dimension |
| `spml.reg` | `spml_reg` | translated projected-normal circular regression likelihood; no Rfast dependency |
| `gcpc.reg` | `gcpc_reg` | translated including fitted anisotropy parameter `rho` |
| `esag.reg` | `esag_reg` | translated 3-D ESAG regression likelihood with joint shape/regression fit |
| `sespc.reg` | `sespc_reg` | translated 3-D SESPC regression likelihood with joint shape/regression fit |
| `spher.reg` | `spher_reg` | translated 3-D spherical-spherical Procrustes rotation regression |
| `embed.aov` | `embed_aov` | arbitrary groups; optional Monte Carlo permutation p-value |
| `hcf.aov` | `hcf_aov` | arbitrary groups; optional correction and Monte Carlo p-value |
| `het.aov` | `het_aov` | arbitrary groups; optional Monte Carlo permutation p-value |
| `lr.aov` | `lr_aov` | common-kappa likelihood-ratio statistic; optional Monte Carlo p-value |
| `embed.boot` | `embed_boot` | centered two-sample bootstrap |
| `hcf.boot` | `hcf_boot` | centered two-sample bootstrap |
| `het.boot` | `het_boot` | centered two-sample bootstrap |
