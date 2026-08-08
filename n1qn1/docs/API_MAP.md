# API map

| R/package concept | Fortran equivalent | Notes |
|---|---|---|
| `n1qn1()` | `n1qn1_minimize()` | Native typed interface |
| `call_eval` | `procedure(n1qn1_objective)` | Scalar objective callback |
| `call_grad` | `procedure(n1qn1_gradient)` | Analytic gradient callback |
| `vars` | `x0(:)` | Initial parameters |
| `epsilon` | `control%epsilon` | Original minimum-coordinate-change tolerance |
| `max_iterations` | `control%max_iterations` | Same role |
| `nsim` | `control%max_evaluations` | Objective/gradient pairs |
| `imp` | `control%verbosity` | Simplified native output |
| `zm`, mode 2 | `initial_hessian(:, :)` | Dense symmetric Hessian, factorized internally |
| `zm`, mode 3 | `initial_factor(:)` | Exact internal LDL-transpose restart |
| `value` | `result%value` | Final objective |
| `par` | `result%x` | Final parameters |
| `H` | `result%hessian` | Reconstructed dense Hessian |
| `c.hess` | `result%c_hess` | Packed dense lower triangle plus zero padding |
| internal `zm(1:nh)` | `result%factor` | Packed LDL-transpose form |
| `n.fn` | `result%function_evaluations` | One per objective/gradient pair |
| `n.gr` | `result%gradient_evaluations` | One per objective/gradient pair |
| `.n1qn1ptr()` | omitted | R dynamic-link interface only |
| `environment`, `assign` | `user_data` and ordinary Fortran assignment | No R environment |
| `print.functions` | callback or `verbosity` | R-specific printing omitted |
