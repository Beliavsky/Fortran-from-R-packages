# Testing

Six deterministic programs cover:

1. generalized inverses, null spaces, empirical multivariate-normal moments,
   KDE, integration, rational approximation, and contrasts;
2. OLS, robust regression, LTS, ridge paths, stepwise AIC, and Box-Cox profiles;
3. distribution fitting, Poisson/NB GLMs, log-linear fitting, and NB simulation;
4. LDA/QDA and robust covariance estimators;
5. correspondence analysis, MCA, classical MDS, Sammon mapping, isoMDS, and
   Shepard disparities;
6. proportional-odds regression, three bandwidth selectors, and
   negative-exponential initialization.

Checked builds enable bounds, allocation, argument, and floating-point runtime
checks through GNU Fortran's `-fcheck=all`.
