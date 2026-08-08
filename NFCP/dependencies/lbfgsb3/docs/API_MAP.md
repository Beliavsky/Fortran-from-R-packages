# R-to-Fortran API map

| R/package component | Fortran replacement | Coverage |
|---|---|---|
| `lbfgsb3c`, `lbfgsb3`, `lbfgsb3f`, `lbfgsb3x` | `lbfgsb_minimize` | Complete computational path with analytic gradient |
| missing `gr`, via `numDeriv::grad` | `lbfgsb_minimize_fd` | Bound-aware finite differences |
| `control$trace` | `lbfgsb_control_t%trace` | Native iteration/evaluation output |
| `control$factr` | `lbfgsb_control_t%factr` | Direct mapping |
| `control$pgtol` | `lbfgsb_control_t%pgtol` | Direct mapping |
| `control$abstol` | `lbfgsb_control_t%abstol` | Direct mapping |
| `control$reltol` | `lbfgsb_control_t%reltol` | Direct mapping |
| `control$lmm` | `lbfgsb_control_t%memory` | Direct mapping |
| `control$maxit` | `max_evaluations` and `max_iterations` | Split into explicit limits |
| `control$info` output | fields of `lbfgsb_result_t` | Always available |
| R `...`/environment | optional polymorphic `user_data` | Native callback state |
| Rcpp function pointers | Fortran procedure callbacks | Native replacement |
| `.lbfgsb3cPtr()` | none | R/C interface-only; omitted |
| `setulb` | `lbfgsb3_core_mod::setulb` | Low-level reverse communication retained |
