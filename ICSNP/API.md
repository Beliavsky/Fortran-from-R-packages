# API reference

All public procedures are available through:

```fortran
use icsnp
```

`dp` is the package real kind.

## Result types

### `type(test_result)`

Fields: `statistic`, `p_value`, `df1`, `df2`, `replications`, `status`, and
`method`.

### `type(location_scatter_result)`

Fields: allocatable `center`, allocatable `scatter`, `iterations`, and
`status`.

### `type(spatial_sign_result)`

Fields: allocatable `signs`, `center`, `shape`, and `status`.

## Pairwise operations

```fortran
call pair_diff(x, result, status)
call pair_sum(x, result, status)
call pair_prod(x, result, status)
```

For an `n x p` matrix, each returns `n*(n-1)/2 x p`, in lexicographic pair
order `(1,2), (1,3), ..., (n-1,n)`.

```fortran
call spatial_ranks(x, ranks)
call signed_ranks(x, ranks)
```

## Location estimators

```fortran
call spatial_median(x, center, status, iterations, init, maxiter, tolerance)
location = hl_loc(x, status)
location = vdw_loc(x, status, int_diff, maxiter)
call HR_Mest(x, result, maxiter, eps_scale, eps_center)
```

## Spatial signs

```fortran
call spatial_sign(x, result, center, shape, estimate_center, estimate_shape, &
   maxiter, tolerance)
```

With neither `center` nor `shape` supplied, both are estimated by the
Hettmansperger-Randles procedure. Logical options permit origin/identity
choices or estimation of only one component.

## Shape and scatter estimators

```fortran
call tyler_shape(x, shape, status, iterations, location, init, maxiter, &
   tolerance, steps)
call duembgen_shape(x, shape, status, iterations, init, maxiter, tolerance, steps)
call duembgen_shape_wt(x, weights, shape, status, iterations, init, maxiter, tolerance)
call symm_huber(x, scatter, status, iterations, qg, init, maxiter, tolerance)
call symm_huber_wt(x, weights, scatter, status, iterations, qg, init, maxiter, tolerance)
call HP1_shape(x, shape, status, location, estimate_location, maxiter, tolerance)
```

Tyler, Duembgen, weighted Duembgen, and HP1 outputs are determinant-normalized.
The Huber procedures return scatter matrices on their natural scale.

## Location tests

```fortran
call HotellingsT2(x, result, y, mu, distribution)
call rank_ctest(x, result, y, mu, scores)
call rank_ctest_groups(x, groups, result, scores)
call rank_ictest(x, result, mu, scores, method, n_simu, seed)
call HP_loc_test(x, result, mu, scores, method, n_perm, seed)
```

Accepted score strings are `sign`, `rank`, and `normal`.

`HotellingsT2` accepts `distribution='f'` or `'chi'`.

`rank_ictest` accepts `method='approximation'`, `'simulation'`, or
`'permutation'`. `HP_loc_test` accepts `approximation` or `permutation`.

## Independence tests

```fortran
call ind_ctest(x, index1, result, index2, scores)
call ind_ictest(x, index1, result, index2, scores, method, n_simu, seed)
```

`index1` and optional `index2` are one-based column-index arrays. If `index2`
is omitted, its columns are the complement of `index1`.

## Status values

- `icsnp_ok`
- `icsnp_invalid_input`
- `icsnp_singular`
- `icsnp_iteration_limit`
- `icsnp_numerical_error`

Use `icsnp_status_message(status)` for a short description.
