# API mapping

| Upstream lmtest API | Fortran computational API | Notes |
|---|---|---|
| `coeftest` | `coefficient_tests` | Takes coefficient vector and covariance matrix directly. |
| `coefci` | `coefficient_confint` | Student-t or normal intervals. |
| `lrtest` | `likelihood_ratio_test` | Pairwise numerical LR kernel; longer sequences are repeated pairwise calls. |
| `waldtest` | `wald_restriction`, `nested_linear_test` | General linear restrictions or nested design matrices. |
| `bgtest` | `breusch_godfrey_test` | LM/F forms; includes auxiliary coefficients and covariance. |
| `bptest` | `breusch_pagan_test` | Studentized/nonstudentized; optional weights and separate variance design `z`. |
| `dwtest` | `durbin_watson_test` | Exact Farebrother AS 153 or upstream normal approximation. |
| `gqtest` | `goldfeld_quandt_test` | Supports point/fraction, alternative, and explicit ordering vector. |
| `harvtest` | `harvey_collier_test` | Recursive-residual implementation. |
| `hmctest` | `harrison_mccabe_test` | Monte Carlo p-value; plotting omitted. `nsim<=0` computes statistic only. |
| `raintest` | `rainbow_test` | Central, explicit-order, and Mahalanobis ordering. |
| `resettest`/`reset` | `reset_test` | `fitted`, `regressor`, and `princomp` variants. |
| `grangertest` | `granger_test` | Matrix/time-series numerical core with lag construction. |
| `coxtest` | `cox_test` | Two design matrices for the two nonnested models. |
| `jtest` | `j_test` | OLS covariance in high-level wrapper; custom covariance can be applied through generic inference routines. |
| `petest` | `pe_test` | Caller supplies the transformed response pair and `islog1/islog2`. |
| `encomptest` | `encompassing_test` | Builds an independent union of the two design spaces and performs nested tests. |
| `pan.f` | `pan_probability` | Free-format module translation; original source retained. |

R-only formula processing, S3 dispatch, model updates, row-name matching, print methods, plotting, and dataset documentation are deliberately outside the numerical port.
