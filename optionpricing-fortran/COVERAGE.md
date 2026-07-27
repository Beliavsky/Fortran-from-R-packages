# Computational coverage

## Exported R routines

| R routine | Fortran routine | Status |
| --- | --- | --- |
| `BS_EC` | `bs_ec`, `bs_european_call` | Implemented; corrected gamma plus legacy field |
| `BS_EP` | `bs_ep`, `bs_european_put` | Implemented; corrected gamma plus legacy field |
| `AsianCall_AppLord` | `asian_call_app_lord` | Implemented |
| `AsianCall` | `asian_call` | Implemented for all working MC/Korobov QMC combinations |

## Analytic and control-variate routines

| R routine | Fortran mapping |
| --- | --- |
| `BS_A` | `bs_a` |
| `Covmat_A` | `covariance_conditional_log_prices` |
| `cond2M_A` | `conditional_average_moments` |
| `evalECV_A` | `eval_ecv_a` |
| `evalECV` | `eval_ecv` |
| `findbcv` | `find_bcv` |
| `evalLB` | `eval_lb` |
| `eval_equad` | `eval_equad` |
| `evalEQCV` | `eval_eqcv` |

## Simulation routines

| R routine | Fortran mapping |
| --- | --- |
| `AsianCall_naive_greeks` | `asian_call_naive_mc` |
| `AsianCall_NCV_LR_greeks` | `asian_call_ncv_lr_mc` |
| `AsianCall_CMC_CV` | `asian_call_cmc_cv` |
| `AsianCall_NCV_CMC_QCV_greeks` | `asian_call_best_mc` |
| `simulateAsianCall_Z` | `simulate_asian_call_z` |
| `AsianCall_naive_greeks_Z` | `asian_call_naive_greeks_z` |
| `AsianCall_NCV_CMC_greeks_Z` | Covered by `conditional_estimates_z` |
| `AsianCall_NCV_CMC_QCV_greeks_Z` | `conditional_estimates_z` and `build_conditional_samples` |

The full conditional/QCV engine returns all six control-variate quantities.
The older price-only and non-QCV variants are implemented by selecting the
needed columns rather than duplicating the same recursion.

## QMC and transformation routines

| R routine | Fortran mapping |
| --- | --- |
| `returnPn_Korobov` | `korobov_lattice` |
| `AsianCall_naive_greeks_qmc` | `asian_call_naive_qmc` |
| `AsianCall_NCV_CMC_greeks_qmc` | Subsumed by `asian_call_best_qmc` |
| `AsianCall_NCV_CMC_QCV_greeks_qmc` | `asian_call_best_qmc` |
| `ortmat` | `orthonormal_complete` internal support |
| `findQmat` | `conditional_generation_matrix` |

The documented generation modes `std`, `pca`, `pcamain`, `lt`, and `ltpca`
are supported for the variance-reduced QMC engine. Naive QMC supports the same
`std` and `pca` modes as the original wrapper.

## Internal numerical infrastructure

- Vector Newton iteration and first-order root initialization
- Bracketed bisection for the Curran boundary
- Adaptive Simpson quadrature for Lord's approximation
- Pivoted QR multivariate regression
- Jacobi symmetric eigendecomposition
- Standard-normal RNG and inverse normal CDF
- Direct, leave-one-out splitting, and pilot-run control-variate regression

## Deliberately excluded

- Histogram and scatter-plot creation
- R list, matrix-name, formula, and `lm` object presentation
- The inactive upstream Sobol branch whose actual Sobol call is commented out

No numerical algorithm depends on R, graphics packages, or external libraries.
