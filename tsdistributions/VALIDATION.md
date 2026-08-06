# Validation

Five deterministic test programs cover:

1. Density/log-density and CDF/quantile consistency for all ten standardized families.
2. Random-generation moments, theoretical moments, and GH/NIG transformations.
3. Numerical maximum-likelihood recovery and Hessian/OPG/sandwich covariance output.
4. Semiparametric GPD-tail/kernel-interior fitting, CDF inversion, density, and RNG.
5. Authorized moment domains, simulation profiling, and invalid-distribution handling.

The build scripts compile with Fortran 2018 conformance, warnings as errors, conversion diagnostics, and full runtime checking in checked mode.
