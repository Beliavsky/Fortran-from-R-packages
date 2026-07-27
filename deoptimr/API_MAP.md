# API map

| Original DEoptimR routine/feature | Modern Fortran implementation | Validation |
|---|---|---|
| `JDEoptim` | `jde_optimize` | Sphere, equality-constrained quadratic, inequality-constrained quadratic, initial-population tests |
| `SPJDEoptim` numerical algorithm | `spjde_optimize` | Aluffi-Pentini, equality and inequality constrained tests |
| `NCDEoptim` | `ncde_optimize` | Four-minimum Becker-Lago, inequality and equality two-minimum tests |
| `handle.bounds` / `handle_bounds` | `handle_bounds` | Exact deterministic unit test |
| Equality conversion `abs(h)-eps` | `transform_constraints` | Scalar/vector behavior tested |
| DE/rand/1/either-or/bin | Internal `reproduce` procedures | Exercised through every optimizer |
| Self-adaptive `F`, `CR`, `pF` | `jde_control`, internal adaptation | Exercised through every optimizer |
| NCDE adaptive neighborhood size | `ncde_control`, internal adaptation | Multimodal tests |
| Median/max stopping reference | `jde_control%compare_to` | Both median and maximum paths tested |
| `add_to_init_pop` | `initial_population` optional matrix | Exact injected-optimum test |
| `details` population output | `de_result` population fields | Size/content tests |
| Constraint archive/population output | `ncde_result` fields | Feasibility and ordering tests |
| Automatic niche radius | negative `ncde_control%niche_radius` | Finite-positive-radius test |
| R `mirai` parallel evaluation | Not translated | External R orchestration, not numerical DE state transition |
| R lists, `...`, argument introspection | Typed callbacks and derived types | Compile-time interfaces |
