# API mapping

Upstream `rmoo` 0.3.2 computational exports map as follows.

| R routine | Fortran counterpart |
|---|---|
| `rmoo`, `nsga` | `rmoo_optimize_*` with `ALG_NSGA1` |
| `nsga2` | `rmoo_optimize_*` with `ALG_NSGA2` |
| `nsga3` | `rmoo_optimize_*` with `ALG_NSGA3` |
| `rnsga2` | `rmoo_optimize_*` with `ALG_RNSGA2` |
| `non_dominated_fronts` | `non_dominated_sort` |
| `crowding_distance` | `crowding_distance` |
| `sharing` | `sharing_dummy_fitness` |
| `modifiedCrowdingDistance` | `rnsga2_survivors` |
| `generate_reference_points` | `generate_reference_points` |
| `get_fixed_rowsum_integer_matrix` | same name |
| `scale_reference_directions` | same name |
| `reference_point_multi_layer` | same name; tensor + row-count API |
| `associate_to_niches` | same name |
| `compute_perpendicular_distance` | same name |
| `compute_niche_count` | same name |
| `niching` | same name |
| `UpdateIdealPoint` | `update_ideal_point` |
| `UpdateWorstPoint` | `update_worst_point` |
| `PerformScalarizing` | `perform_scalarizing` |
| `get_nadir_point` | same name |
| `calc_norm_pref_distance` | same name |
| `generational_distance` | same name |
| GD+/IGD utilities | `gd_plus`, `igd`, `igd_plus` |
| real population | vendored GA `random_real_population` / high-level initializer |
| binary population | vendored GA `random_binary_population` / high-level initializer |
| permutation population | vendored GA `random_perm_population` / high-level initializer |
| discrete population | high-level integer initializer |
| `rmoo_tourSelection` | `tournament_nsga1`, `tournament_nsga2`, `tournament_rank` |
| `rmoo_lrSelection` | `linear_rank_selection` |
| `rmooreal_sbxCrossover` | `sbx_crossover` |
| `rmoo_spCrossover` | `single_point_crossover_real/int` |
| `rmoo_uxCrossover` | `uniform_crossover_real/int` |
| `rmoo_huxCrossover` | `hux_crossover` |
| `rmooperm_oxCrossover` | `ox_crossover` |
| `rmooreal_polMutation` | `polynomial_mutation` |
| `rmooreal_raMutation` | `random_real_mutation` |
| `rmoobin_raMutation` | `random_binary_mutation` |
| `rmooperm_simMutation` | `inversion_mutation` |
| `rmoo_uxMutation` | `uniform_integer_mutation` |

R object getters (`getPopulation`, `getFitness`, `getCrowdingDistance`, etc.)
are unnecessary because the corresponding arrays are public components of
`rmoo_real_result` and `rmoo_integer_result`.

`rmooControl` is represented by optional high-level arguments.  Arbitrary R
operator callbacks are not mirrored in v0.1.0; the translated low-level
operators are public so a caller can assemble a custom evolutionary loop.

`startParallel`, `stopParallel`, plotting, monitor, print, progress and summary
methods are R/UI infrastructure and are intentionally omitted.
