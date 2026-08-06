# Translation coverage

## Directly represented computational behavior

- Bisquare and Huber rho/psi behavior, numerical second derivatives, robust weights, Gaussian-efficiency calculation, tuning search, M-scale iteration, and inverse robust R-squared.
- DCML coefficient mixing, scale recomputation, covariance construction, RFPE decomposition, and robust nested-model test statistic.
- BY logistic fitting through the supplied robustbase translation; WBY and WML leverage screening and full-data predictions.
- Classical covariance, MM-SHR and Rocke descent, robust distances, correlation conversion, and consistency rescaling.
- Residual M-scale PCA iteration: Locantore initialization, weighted center/scatter updates, eigensystem updates, residual scales, loading convergence, scores, and fitted values.

## Self-contained numerical equivalents

- The R package obtains Pen~a-Yohai candidates from `pyinit`. This port uses the supplied robustbase S-regression implementation as the high-breakdown start, followed by the RobStatTM fixed-scale M refinement.
- `SMPY` in R separates factor and continuous columns using model-frame metadata. The Fortran matrix API has no factor metadata; `sm_py_fit` uses the same robust numerical path on the supplied design matrix.
- `KurtSDNew` contains extensive directional optimization code. `init_pp` retains its computational role using deterministic data directions plus seeded random unit directions and robust projection outlyingness.
- `fastmve` is supplied by the vendored rrcov Fortran implementation rather than by translating the package's C entry point line-for-line.
- `step.lmrobdetMM` becomes `stepwise_rfpe`, operating on matrix columns rather than R formula terms and hierarchy attributes.
- `optv0` and `moptv0` expose scalar tuning constants in Fortran rather than R's named tuning vectors. The scalar loss API uses the optimal redescending family represented by the numerical kernels.

## Omitted non-numerical infrastructure

- Plotting and graphics.
- R formulas, terms objects, model frames, contrasts, factor-level metadata, and `na.action` dispatch.
- S3/S4 constructors, methods, `print`, `summary`, `predict`, and class attributes.
- R package datasets as callable Fortran data modules.
- R registration, `.Call`, `.Fortran`, and R memory-management wrappers.
- Documentation rendering and examples requiring the R runtime.

The unmodified upstream package snapshot is retained in `upstream/RobStatTM-master`.
