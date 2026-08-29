# API mapping

| Upstream computational area | Fortran API |
|---|---|
| `phdensity`, `phcdf`, PH `dens/cdf/haz/quan/moment/mean/var/laplace/mgf` | `matrixdist_ph`: `ph_density`, `ph_cdf`, `ph_survival`, `ph_hazard`, `ph_quantile`, `ph_moment`, `ph_mean`, `ph_variance`, `ph_laplace`, `ph_mgf` |
| PH `sim`, `logLik`, EM | `rphasetype`, `ph_loglik`, `emstep_ph`, `fit_ph_em` |
| `sum_ph`, minimum/maximum/mixture | `ph_sum`, `ph_minimum`, `ph_maximum`, `ph_mixture` |
| `mweibull*`, `mpareto*`, `mlognormal*`, `mloglogistic*`, `mgompertz*`, `mgev*`, `riph`, `rmatrixgev` | `matrixdist_iph` and `matrixdist_simulation` |
| `dphdensity`, `dphcdf`, `dph_pgf`, moments/simulation/EM | `matrixdist_dph`, `rdphasetype`, `emstep_dph`, `fit_dph_em` |
| `bivph_density`, `bivph_tail`, `bivph_laplace`, moments/simulation/EM | `matrixdist_multivariate`, `rbivph`, `matrixdist_bivariate_fit` |
| `bivdph_density`, `bivdph_tail`, PGF/moments/simulation/EM | `matrixdist_multivariate`, `rbivdph`, `matrixdist_bivariate_fit` |
| `mdphdensity`, mDPH PGF/moments/simulation/EM | `mdph_density_point`, `mdph_pgf_point`, `mdph_factorial_moment`, `rmdph`, `emstep_mdph` |
| mPH density/CDF/Laplace/moments/right-censored EM | `mph_density_point`, `mph_cdf_point`, `mph_laplace_point`, `mph_moment`, `mph_loglik`, `rmph`, `emstep_mph_rc` |
| `biviph`, `miph` fixed transforms | `biviph_density`, `biviph_tail`, `miph_density_point`, `miph_cdf_point`, `rbiviph`, `rmiph` |
| MPH* mean/variance/simulation and reward transform | `mphstar_mean`, `mphstar_cov`, `rmphstar`, `tvr_ph`, `linear_combination` |
| `random_reward`, `rew_sanity_check`, `marginal_expectation`, MPH* EM | `random_reward`, `reward_sanity_check`, `marginal_sojourn_expectation`, `emstep_mphstar` |
| `tvr_ph`, `tvr_dph`, `merge_matrices`, `plus_states` | `matrixdist_transformations` |
| `random_structure`, `random_structure_bivph` | `matrixdist_structure` |
| `matrix_exponential`, `matrix_vanloan`, `matrix_power`, `kronecker_sum`, inverse/product helpers | `matrixdist_linalg` |
| `find_n`, `pow2_matrix`, uniformization-related helpers | `matrixdist_numerics` |
| survival PH `evaluate`/fixed-parameter likelihood | `sph_density`, `sph_survival`, `sph_loglik` |
| `LRT`, covariance-to-correlation behavior | `lrt`, `cov_to_cor` |

### Consolidated/omitted interfaces

The upstream RK, Pade and uniformization versions of the same EM/likelihood calculation are not duplicated as separate public entry points.  The port retains exact Pade/Van-Loan kernels plus uniformization helpers.  R formula parsing, S4 constructors/accessors, plotting/printing, profile/progress output, and `nnet::multinom` formula orchestration are omitted.
