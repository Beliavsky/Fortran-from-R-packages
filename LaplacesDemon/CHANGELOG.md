# Changelog

## 0.3.1

- Fixed `test_v02` portability across Fortran runtimes with different intrinsic
  `random_number` streams.
- Replaced brittle NUTS and pCN Monte Carlo thresholds with statistically
  symmetric/robust checks; the sampler implementations are unchanged.
- Added named assertion diagnostics to `test_v02`, so future failures identify
  the exact sampler/statistic instead of only reporting `FAIL 1`.
- Stress-tested the revised v0.2 sampler assertions over 50 alternate seeds
  with zero failures.

## 0.3.0

- Completed the fixed-dimensional LaplacesDemon sampler catalog with ADMG,
  AFSS, adaptive/griddy Gibbs, AHMC, AMM/AMWG, CHARM/HARM, DRM, IM, INCA,
  MCMCMC, MTM, OHSS, RDMH, refractive sampling, RSS, SAMWG/SMWG and updating
  variants, SGLD, THMC, UESS, callback Gibbs and fixed-vector RJ variable
  selection.
- Corrected RDMH to the upstream `epsilon = U^s`, `s in {-1,+1}` proposal and
  Jacobian rather than a lognormal multiplicative approximation.
- Added the remaining `LaplaceApproximation` optimizer catalog: AGA, BHHH, CG,
  DFP, HAR, HJ, LBFGS, LM, NM, NR, PSO, Rprop, SGD, SOMA, SPG, SR1 and TR.
- Added componentwise iterative quadrature and adaptive sparse-grid quadrature.
- Added BMK, Heidelberger-Welch and KS diagnostics.
- Expanded the distribution/prior catalog with matrix gamma, inverse matrix
  gamma, multivariate Cauchy/Laplace/power-exponential variants, mixture and
  shrinkage priors, Wishart/inverse-Wishart RNG/Cholesky APIs, generic
  truncation and scalar compatibility CDF/quantile/RNG routines.
- Added the remaining exported `distributions.R` kernels: CRMRF,
  Huang-Wand-Cholesky and Laplace-mixture density alias.
- Added `test_v03`, `test_v03_completion`, `test_distributions_v03` and
  `example/completion_v03.f90`.
- Added selected upstream v0.3 R source files for provenance.

## 0.2.0

- Added NUTS with recursive tree doubling and dual-averaged step-size adaptation.
- Added HMC dual averaging (HMCDA).
- Added DRAM, RAM, pCN and the two-state t-walk.
- Added the Salimans2 variational-Bayes stochastic approximation.
- Added Raftery-Lewis and Hangartner diagnostics.
- Added `test_v02` and expanded API/translation documentation.

## 0.1.0

- Initial modern Fortran/FPM computational port of LaplacesDemon 16.1.8.
- Added numerical differentiation, SPD linear algebra, BFGS and Laplace
  approximation.
- Added tensor adaptive Gauss-Hermite approximation and quadrature rules.
- Added RWM, MWG, Adaptive Metropolis, MALA, HMC, slice, elliptical slice,
  AIES and DEMC samplers.
- Added normal importance sampling, SIR, rejection sampling, PMC and Bayesian
  bootstrap weights.
- Added IAT/ESS/MCSE/Geweke/Gelman/WAIC/KLD diagnostics and marginal-likelihood
  helpers.
- Added representative scalar/multivariate distributions and priors.
- Preserved upstream MIT licensing/provenance and selected reference R files.
