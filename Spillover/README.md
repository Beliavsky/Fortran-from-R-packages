# Spillover modern Fortran

This package translates the computational core of the R package **Spillover 0.1.1** to modern Fortran 2018. It provides self-contained VAR estimation, forecast-error variance decompositions, Diebold-Yilmaz connectedness tables, and rolling connectedness calculations. It uses no external numerical library.

## Implemented functionality

- OLS estimation of VAR(p) models with no deterministic term, a constant, a trend, or both
- VAR moving-average coefficient recursion corresponding to `vars::Phi`
- Generalized FEVD of Pesaran and Shin, normalized or unnormalized
- Cholesky-orthogonalized FEVD
- Single-order, sampled-permutation, and exact-permutation orthogonalized tables
- Generalized and orthogonalized total, directional, and net connectedness
- Compatibility tables corresponding to `G.spillover()` and `O.spillover()`
- Rolling total and net connectedness
- Full rolling generalized directional connectedness, including pairwise net measures
- Compatibility wrappers named `g_fevd`, `g_spillover`, `o_spillover`, `roll_spillover`, `roll_net`, and `total_dynamic_spillover`

Plotting, `zoo` indexing, R data frames, S3 objects, and R package infrastructure are intentionally omitted.

## Build

With FPM:

```text
fpm test
fpm run --example example_spillover
```

With GNU Make:

```text
make check
make release
```

`make check` uses strict Fortran 2018 compilation, warnings as errors, bounds and runtime checking, and backtraces. `make release` uses `-O3` while retaining warnings as errors.

## Minimal use

```fortran
use spillover

type(var_model) :: model
type(spillover_result) :: connectedness
integer :: info
character(len=200) :: message

call fit_var(returns, 2, model, var_const, info, message)
call g_spillover(model, 10, .true., connectedness, info, message)
print *, connectedness%total
```

The matrix `connectedness%shares(i,j)` is the percentage contribution from shock `j` to the forecast-error variance of variable `i`. Therefore:

- `from(i)` is the off-diagonal row sum for variable `i`;
- `to(i)` is the off-diagonal column sum;
- `net(i) = to(i) - from(i)`;
- `total` is the sum of all off-diagonal shares divided by the number of variables.

## Orthogonalized ordering averages

`ortho_single` uses the original variable order. `ortho_partial` samples permutations, and `ortho_total` enumerates all permutations. Exact enumeration is intentionally guarded by a configurable dimension limit, because the number of permutations is factorial. The original `O.spillover()` source accidentally interchanges the `partial` and `total` calls; pass `source_compatible=.true.` to reproduce that behavior.

## License

The upstream package declares `GPL-2`. This translation is distributed under **GPL-2.0-only**. The original source archive is retained under `upstream/`, and the complete GPL version 2 text is included in `LICENSE` and `license/GPL-2.0.txt`.
