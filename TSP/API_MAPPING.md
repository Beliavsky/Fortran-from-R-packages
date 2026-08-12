# R-to-Fortran computational API mapping

| Upstream R/C functionality | Fortran |
|---|---|
| `TSP`, `ATSP` numeric representation | square `real(dp)` cost matrix |
| `ETSP`, `as.TSP.ETSP` | `euclidean_distance_matrix` |
| `TOUR` | `type(tsp_tour)` or integer permutation |
| `tour_length` / C kernels | `tour_length`, `etsp_tour_length` |
| C `insertion_cost` | `insertion_cost` |
| `tsp_nn` | `nearest_neighbor` |
| `tsp_repetitive_nn` | `repetitive_nearest_neighbor` |
| `tsp_insertion` | `insertion_heuristic` |
| `tsp_insertion_arbitrary` | `arbitrary_insertion` |
| C/R `two_opt` | `two_opt`; registered symmetric kernel: `two_opt_symmetric` |
| `tsp_SA` and local moves | `simulated_annealing`, `sa_reversal`, `sa_swap`, `sa_mixed` |
| `solve_TSP` | `solve_tsp` |
| `.replaceInf` | `replace_infinite` |
| `insert_dummy` | `insert_dummy` |
| `reformulate_ATSP_as_TSP` | `reformulate_atsp_as_tsp` |
| `filter_ATSP_as_TSP_dummies` | `filter_atsp_tour` |
| `cut_tour` | `cut_tour_single`, `cut_tour_multiple` |
| `read_TSPLIB` | `read_tsplib` |
| `write_TSPLIB.*` | `write_tsplib_tsp`, `write_tsplib_atsp`, `write_tsplib_etsp` |
| `.tsplib_att_dist` | `tsplib_att_distance` |
| `.tsplib_geo_dist` | `tsplib_geo_distance` |
| plotting/image methods | omitted |
| Concorde/linkern interfaces | omitted; external solver code is not upstream |

Solver constants such as `tsp_nn`, `tsp_cheapest_insertion`, and
`tsp_sa_method` are defined in `tsp_types` and re-exported by the umbrella
`tsp` module. Controls are collected in `type(tsp_control)`.
