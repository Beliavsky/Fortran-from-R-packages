# Translation coverage

| Upstream area | Fortran coverage |
|---|---|
| `maxLik` dispatcher | Common `max_lik` method dispatcher |
| `maxNR` / `maxNRCompute` | Newton-Raphson, step halving, Hessian correction, numeric derivatives |
| `maxBFGS` | BFGS maximization |
| `maxBFGSR` / `maxBFGSRCompute` | BFGSR-compatible entry point using the self-contained BFGS engine |
| `maxBHHH` | Observation-score BHHH/Fisher scoring |
| `maxCG` | Polak-Ribiere nonlinear conjugate gradient |
| `maxNM` | Nelder-Mead simplex maximization |
| `maxSANN` | Seeded simulated annealing |
| `maxSGA` / `maxSGACompute` | Mini-batch stochastic-gradient ascent, momentum, clipping, patience |
| `maxAdam` | Mini-batch Adam with bias correction |
| `numericGradient` | Scalar objective numerical gradient plus general vector numerical Jacobian |
| `numericHessian` / `numericNHessian` | Numerical Hessian from analytic or numerical gradients |
| `compareDerivatives` | Analytic/numerical gradient and Hessian reports |
| `sumt` | Outer quadratic-penalty sequence for linear equality constraints |
| `constrOptim2` | Linear inequality constraints through the same penalty framework |
| `activePar`, `prepareFixed`, `addFixedPar`, `fnSubset` | Active mask, fixed indices, packing and expansion utilities |
| `vcov`, `stdEr`, `bread`, `estfun` | Hessian covariance, robust score covariance, standard errors, stored scores |
| `condiNumber` | Hessian and progressive matrix condition numbers |
| `AIC`, `confint` | AIC and normal confidence-interval utilities |
| stored values/parameters | Optional objective and parameter histories |
| R S3/S4 classes and methods | Replaced by Fortran derived types |
| `summary`, `print`, `tidy`, `glance` | Not translated; presentation infrastructure |
| Formula/model-frame and `...` handling | Not applicable to explicit Fortran callbacks |
| Vignettes and package build machinery | Retained upstream, not translated |
| Plotting | No plotting code in scope |

The complete upstream R package used for the port is retained under
`original/maxLik-master/`.
