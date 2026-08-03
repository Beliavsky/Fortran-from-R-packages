# Translation coverage

The original `NAMESPACE` exports 24 computational names. All 24 are represented
in the public Fortran API.

| Original export | Fortran API |
|---|---|
| `A`, `Astar`, `Ainv` support | `A`, `Astar`, `Ainv` |
| `D`, `Dstar` | `D`, `Dstar` |
| `L`, `Lstar`, `Linv` support | `L`, `Lstar`, `Linv` |
| `accuracy` | `accuracy` |
| `block_diag` | generic `block_diag` for two or three matrices |
| `cluster_k_component_graph` | same name, subroutine |
| `fdr`, `fscore`, `npv`, `recall`, `specificity` | same names |
| `learn_bipartite_graph` | same name, subroutine |
| `learn_bipartite_k_component_graph` | same name, subroutine |
| `learn_combinatorial_graph_laplacian` | same name, subroutine |
| `learn_graph_sigrep` | same name, subroutine |
| `learn_k_component_graph` | same name, subroutine |
| `learn_laplacian_gle_admm` | same name, subroutine |
| `learn_laplacian_gle_mm` | same name, subroutine |
| `learn_smooth_approx_graph` | same name, subroutine |
| `learn_smooth_graph` | same name, subroutine |
| `relative_error` | same name |

Additional translated computational routines include `learn_cospectral_graph`,
`Mmat`, `Pmat`, `Dmat`, `vec`, `vecLmat`, `pairwise_matrix_rownorm2`,
`upper_view_vec`, `metrics`, and `prial`.

Plotting shown in R examples and vignettes is external to the original package's
computational implementation and is omitted. See `original/OMITTED.md`.
