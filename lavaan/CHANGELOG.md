# Changelog

## 0.7.0

- Added `muthen1984_mixed` for mixed continuous/ordinal statistics with optional conditional-x regressions.
- Added numeric-numeric, numeric-ordinal (polyserial), and ordinal-ordinal conditional latent-correlation estimation.
- Added optional delete-one stagewise Gamma/WLS construction for the mixed Muthen statistic vector.
- Added missing-response random-coefficient likelihood, ML fitting, and empirical-Bayes random-effect prediction.
- Added automatic marker selection and name-aware MIIV equation output, plus a small parameter-table marker descriptor.
- Added SAM joint stage-1/stage-2 block covariance construction and generic second-order delta bias propagation.
- Added explicit Satorra-Bentler 2001 and 2010 scaled chi-square difference-test calculations.
- Added v0.7 parity tests for all new paths and preserved all v0.1-v0.6 tests.

## 0.6.0

- Added the all-ordinal Muthen (1984) three-stage NACOV path with analytic threshold/polychoric scores and explicit A11/A21/A22 block bread.
- Added `categorical_wls_statistics_muthen` and detailed Muthen diagnostic result matrices.
- Added marker-variable MIIV rewriting with loading-ratio coefficient conversion and composite-disturbance instrument screening.
- Added posterior-adaptive Gauss-Hermite mixed-outcome MML and fitting for one to three latent factors.
- Added continuous ADF SAM Gamma and Browne finite-sample unbiased Gamma calculations.
- Added general Gaussian random-coefficient/random-slope ML with unrestricted G/R covariance matrices and empirical-Bayes random-effect prediction.
- Added Browne residual normal-theory/ADF test utilities.
- Added v0.6 parity tests with closed-form binary Muthen NACOV, adaptive-Gaussian MML, MIIV marker rewriting, random-effect likelihood/EB, SAM Gamma, and Browne-test checks.

## 0.5.0

- Added exact Hayakawa corrected `tr(U Gamma)` / `tr((U Gamma)^2)` estimation and corrected adjusted-test formulas.
- Added ULS, GLS, 2RLS and RLS covariance-statistic MIIV estimators and their upstream-style Jacobian linearizations.
- Added `miiv_vcov_from_gamma` for asymptotic covariance propagation.
- Added matrix-level Yuan-Chan global SAM scaling through residualized measurement-statistic covariance.
- Added deterministic Halton-QMC mixed continuous/ordinal MML for up to eight latent factors.
- Added v0.5 parity tests including finite-difference MIIV Jacobian checks and a four-factor known-marginal MML check.

## 0.4.0

- Added local SAM stage-1 uncertainty propagation with delta-method stage-2 covariance correction.
- Added automatic RAM MIIV equation discovery, descendant exclusion, relevance screening, and identification flags.
- Added influence-function threshold/polychoric categorical Gamma and WLS/DWLS weights without jackknife refits.
- Added exact missing-data two-level Gaussian FIML and unrestricted within/between H1 likelihood-ratio testing.
- Added mixed continuous/ordinal Gauss-Hermite MML with missing responses and arbitrary fixed latent mean/covariance.
- Added scaled-shifted and Yuan-Bentler robust-test calculations.
- Added v0.4 parity tests for all refinement areas.

## 0.3.0

- Added principal-axis and Gaussian-ML EFA extraction.
- Vendored GPArotation-fortran and added orthogonal/oblique EFA rotation.
- Added SAM-style two-stage covariance/raw-data fitting and measurement-parameter freezing.
- Added raw-data and covariance-form MIIV/2SLS, first-stage F diagnostics, Sargan statistic, and RAM MIIV candidate screening.
- Added threshold-inclusive categorical statistic vectors and jackknife Gamma/WLS/DWLS matrices.
- Added exact complete-data two-level Gaussian RAM likelihood fitting.
- Added Gauss-Hermite ordinal factor marginal maximum likelihood for one to three latent dimensions.
- Added mixed continuous/ordinal pairwise likelihood.
- Added mean and mean/variance adjusted covariance-structure robust test corrections.
- Added v0.3 parity tests covering the new workflows.

## 0.2.0

- Added linked-parameter multigroup covariance/mean estimation.
- Added raw-data multigroup FIML with missing-data support.
- Added nonlinear equality/inequality constrained covariance fitting and tangent-space vcov correction.
- Added casewise and cluster-robust sandwich inference.
- Added covariance-structure Satorra-Bentler scaling.
- Added nuisance-adjusted modification indices, EPCs, and joint score tests.
- Added automatic jackknife WLS/DWLS weights for polychoric correlation structures.
- Added ordinal pairwise-likelihood RAM fitting.
- Added nonparametric bootstrap SE/bias/percentile intervals.
- Added within/between clustered moment decomposition, separate RAM level fits, and ICCs.

## 0.1.0

- Initial RAM/LISREL SEM core with ML/GLS/ULS/WLS/DWLS, FIML, ordinal polychorics, simulation, fit diagnostics, and factor scores.
