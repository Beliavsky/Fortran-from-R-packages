# Validation

The project is validated with direct GNU Fortran builds using Fortran 2018, runtime checking, and `-Werror=implicit-interface`, linked against BLAS/LAPACK.

Retained tests cover:

- one-state exponential PH identities for density/CDF/moments/Laplace/quantiles;
- two-state PH values independently computed with SciPy `expm`/linear algebra;
- one- and two-state DPH identities;
- all six inhomogeneous transforms plus exact one-state GEV consistency;
- feed-forward bivariate PH/DPH densities, means and covariance;
- conditional multivariate PH/DPH densities and moments;
- MPH* moments/reward covariance;
- Pade versus uniformization matrix-exponential agreement;
- PH, DPH, bivariate PH, bivariate DPH, mPH, mDPH and MPH* EM one-state MLE identities;
- arbitrary right-censoring in mPH (`events / total exposure` in the one-state case);
- PH/DPH and bivariate simulation sample means;
- survival-regression `reg` and `aft` special cases;
- reward time scaling and reward-row normalization.

FPM itself was not available in the execution environment, so the FPM source graph is checked with an equivalent dependency-ordered direct compiler build.
