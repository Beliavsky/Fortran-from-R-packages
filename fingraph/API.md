# API

All public numerical definitions are available from:

```fortran
use fingraph
```

For narrower imports, use `only` clauses.

## Kinds, status, and result type

- `dp = kind(1.0d0)`
- `type(fingraph_result)`
- `fg_ok`
- `fg_invalid_input`
- `fg_size_mismatch`
- `fg_singular_matrix`
- `fg_no_convergence`
- `fg_not_positive_definite`
- `fg_allocation_failure`
- `fg_infeasible`
- `fg_status_message(status)`

`fingraph_result` contains:

- `laplacian`, `adjacency`, `theta`, and `weights`
- `primal_lap_residual`, `primal_deg_residual`, and `dual_residual`
- `lagrangian`, `elapsed_time`, and `beta_seq`
- `convergence`, `iterations`, `status`, `rho`, and `beta`

## Graph estimators

### `learn_connected_graph`

```fortran
call learn_connected_graph(s,result,initial_weights,initialization, &
   degree,degrees,rho,maxiter,reltol)
```

Required arguments:

- `s(:,:)`: square covariance matrix.
- `result`: `type(fingraph_result)` output.

Defaults:

- `initialization='naive'`
- node degree `1`
- `rho=1`
- `maxiter=10000`
- `reltol=1e-5`

`degrees(:)` supplies one target degree per node. `degree` supplies one common
target degree. Do not pass both. `initial_weights(:)` overrides the named
initialization method.

### `learn_regular_heavytail_graph`

```fortran
call learn_regular_heavytail_graph(x,result,heavy_type,nu, &
   initial_weights,initialization,degree,degrees,rho,update_rho, &
   maxiter,reltol)
```

- `x(n,p)` contains observations by rows and graph nodes by columns.
- `heavy_type` is `'gaussian'` or `'student'`.
- Student-t fitting requires `nu > 2`.
- `update_rho` defaults to `.true.`.

The output includes residual, augmented-Lagrangian, and elapsed-time histories.

### `learn_kcomp_heavytail_graph`

```fortran
call learn_kcomp_heavytail_graph(x,result,k,heavy_type,nu, &
   initial_weights,initialization,degree,degrees,beta,update_beta, &
   early_stopping,rho,update_rho,maxiter,reltol,record_objective)
```

Defaults mirror the R function:

- `k=1`
- `heavy_type='gaussian'`
- `beta=1e-8`
- `update_beta=.true.`
- `early_stopping=.false.`
- `rho=1`
- `update_rho=.false.`
- `maxiter=10000`
- `reltol=1e-5`
- `record_objective=.false.`

The observations are centered and standardized internally, matching R's
`scale(as.matrix(X))` call.

## Internal computational helpers exposed for testing

```fortran
weight = compute_student_weight(w,lstar_sq,p,nu)
value = compute_augmented_lagrangian_ht(...)
value = compute_augmented_lagrangian_kcomp_ht(...)
```

The two augmented-Lagrangian routines use a two-dimensional `lstar_sq(n,m)`
array, where row `q` stores the edge-coordinate representation of observation
`q`'s outer product.

## Reused graph API

The following graph routines are also public:

- `L`, `A`, `D`
- `Lstar`, `Astar`, `Dstar`
- `Linv`, `Ainv`
- `Mmat`, `Pmat`, `Dmat`, `vec`, `vecLmat`
- `relative_error`, `fscore`, `metrics`, `block_diag`

## Simulation support

The separate `fingraph_rng` module provides deterministic portable test and
example generators:

```fortran
use fingraph_rng, only : random_mvn, random_mvt
```

These are not part of the original Fingraph export list.
