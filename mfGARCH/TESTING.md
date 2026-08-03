# Testing

The regression suite covers:

1. beta-MIDAS weights, long-run variance construction, GJR recursion, and
   parameter transformations;
2. likelihood identities and multi-horizon prediction;
3. reproducible Gaussian and Student-t simulation;
4. simple-GARCH estimation and all three covariance estimators;
5. GARCH-MIDAS estimation with a low-frequency covariate;
6. two-covariate long-run components;
7. realized-variance-dependent and diffusion-limit simulation.

The supplied scripts compile every module with warnings enabled, run all test
programs, and compile/run all examples and the demonstration program. The
strict script additionally enables runtime bounds, allocation, and argument
checks.
