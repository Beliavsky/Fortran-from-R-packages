# API map

| roptim C++ / R concept | Fortran API | Notes |
|---|---|---|
| `Roptim<Functor>` | `roptim_minimize` | Unified procedure with callbacks |
| `Functor::operator()` | objective procedure | Required |
| `Functor::Gradient` | optional gradient procedure | Finite differences if absent |
| `Functor::ApproximateGradient` | `roptim_approximate_gradient` | Central, bound-aware differences |
| `Functor::ApproximateHessian` | `roptim_approximate_hessian` | Symmetrized gradient differences |
| `RoptimControl` | `roptim_control_t` | Typed control values |
| `par`, `value`, counts, convergence | `roptim_result_t` | Allocatable parameter and Hessian fields |
| `set_lower`, `set_upper` | `lower=`, `upper=` | Scalar or elementwise arrays |
| `set_hessian(true)` | `control%compute_hessian=.true.` | Computed after optimization |
| `Nelder-Mead` | `method_nelder_mead` | Native implementation |
| `BFGS` | `method_bfgs` | Full inverse-BFGS update |
| `CG` | `method_cg` | `cg_type=1,2,3` |
| `L-BFGS-B` | `method_lbfgsb` | L-BFGS-B 3.0 reverse communication |
| `SANN` | `method_sann` | Logarithmic cooling schedule |
| custom SANN generator | optional `proposal=` callback | Explicit candidate callback |
| C++/R printing | ordinary Fortran I/O | No class-printing subsystem |
