# API overview

All public entities are available through:

```fortran
use spectral_graph_topology
```

## Kinds, status, and results

- `dp`: double-precision real kind.
- `graph_result`: matrices, weights, spectral variables, histories,
  convergence flag, iteration count, and status.
- `graph_metrics`: F-score, recall, specificity, accuracy, NPV, FDR, and counts.
- `sgt_ok`, `sgt_invalid_input`, `sgt_size_mismatch`,
  `sgt_singular_matrix`, `sgt_no_convergence`,
  `sgt_not_positive_definite`, and `sgt_infeasible`.
- `sgt_status_message(status)`.

## Operators

The edge vector ordering is `(1,2), (1,3), ..., (1,n), (2,3), ...`.

- `L(w)`: edge weights to a combinatorial Laplacian.
- `A(w)`: edge weights to a symmetric adjacency matrix.
- `D(w)`: weighted node degrees.
- `Lstar(M)`, `Astar(M)`, `Dstar(v)`: adjoint operators.
- `Linv(M)`, `Ainv(M)`: matrix-to-edge-vector inverses.
- `Mmat(m)`, `Pmat(m)`, `Dmat(m)`: matrix representations of operator
  compositions, where `m` is an edge-vector length.
- `vec(M)`, `vecLmat(n)`.

## Utilities and metrics

- `block_diag(A,B)` and `block_diag(A,B,C)`.
- `relative_error(West,Wtrue)`.
- `metrics(Wtrue,West,eps)`.
- `fscore`, `recall`, `specificity`, `fdr`, `npv`, `accuracy`.
- `pairwise_matrix_rownorm2(M)`.
- `upper_view_vec(M)`.
- `prial(Ltrue,Lest,Lscm)`.

## Learning routines

Each routine has the form:

```fortran
call routine_name(input, result, optional_arguments...)
```

where `result` is `type(graph_result)`.

### Spectral constraints

- `learn_k_component_graph(S,result,...)`
- `learn_cospectral_graph(S,fixed_lambda,result,...)`
- `learn_bipartite_graph(S,result,...)`
- `learn_bipartite_k_component_graph(S,result,...)`

For the routines accepting either covariance or data, set
`is_data_matrix=.true.` for a `p x observations/features` matrix. A nonsquare
input is automatically treated as data. Set `use_qp=.true.` for nonnegative
least-squares initialization.

### Smooth graphs and clustering

- `learn_smooth_approx_graph(Y,m,result)`
- `cluster_k_component_graph(Y,result,...)`
- `learn_smooth_graph(X,result,...)`
- `learn_graph_sigrep(X,result,...)`

### Laplacian estimators

- `learn_laplacian_gle_mm(S,result,...)`
- `learn_laplacian_gle_admm(S,result,...)`
- `learn_combinatorial_graph_laplacian(S,result,...)`

The optional `a_mask` is a square real matrix whose positive off-diagonal entries
mark allowed edges.

## `graph_result` fields

Fields are allocated only when relevant:

- `laplacian`, `adjacency`, `weights`
- `eigenvalues`, `eigenvectors`
- `auxiliary_eigenvalues`, `auxiliary_eigenvectors`
- `smoothed_data`
- `objective`, `negative_log_likelihood`
- `elapsed_time`, `parameter_history`, `weight_history`
- `convergence`, `iterations`, `status`
- final `beta`, `nu`, and `lipschitz` values
