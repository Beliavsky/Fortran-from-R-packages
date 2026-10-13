# API and coverage guide

Import from `entropy_module`. Every real argument uses the public `dp` kind.
Arguments in square brackets below are optional; they can be passed by keyword.
All numerical functions/subroutines are pure. The README's full 33-entry table
and `fpm.toml` form the matching exact R-to-Fortran coverage inventory.

## Frequencies and entropy

- `freqs_empirical(y)` normalizes a vector or matrix.
- `freqs_dirichlet(y, a)` normalizes counts plus pseudocounts. `y` is a vector
  or matrix; `a` is a scalar or an array with exactly the same shape.
- `freqs_shrink(y [, lambda])` blends empirical frequencies with a uniform
  target over every element. It preserves vector/matrix shape.
- `shrinkage_intensity(y)` gives the estimated lambda for a count vector.
  This useful public helper does not add an R coverage entry because its
  upstream counterpart is internal.
- `freqs(y [, method, lambda])` chooses an estimator for a count vector.
- `entropy_plugin(weights [, unit])` computes normalized plug-in entropy.
- `entropy_empirical(y [, unit])` estimates entropy from a count vector.
- `entropy_dirichlet(y, a [, unit])` uses a scalar per-bin pseudocount.
- `entropy_shrink(y [, lambda, unit])` estimates shrinkage entropy.
- `entropy_millermadow(y [, unit])` adds the observed-bin bias correction.
- `entropy_chaoshen(y [, unit])` applies coverage and detection correction.
- `entropy(y [, method, unit, lambda])` dispatches among the eight implemented
  methods. Full case-sensitive method strings are listed in README.md.

For a matrix entropy, explicitly flatten using `reshape(y, [size(y)])`.
For cell-specific Dirichlet priors, use
`entropy_plugin(freqs_dirichlet(y, a))` with vector inputs.

## Divergences

`op` below means either `kl` or `chi2`:

- `op_plugin(weights1, weights2 [, unit])`
- `op_empirical(y1, y2 [, unit])`
- `op_dirichlet(y1, y2, a1, a2 [, unit])`
- `op_shrink(y1, y2 [, lambda1, lambda2, unit])`

Both sample vectors must have equal lengths. Their totals need not match;
each is separately normalized. Optional lambdas are independently estimated
when absent. Dirichlet wrapper priors are scalar per-bin pseudocounts.

## Two-variable information and independence

`op` below means either `mi` or `chi2indep`:

- `op_plugin(weights [, unit])`
- `op_empirical(y [, unit])`
- `op_dirichlet(y, a [, unit])`
- `op_shrink(y [, lambda, unit])`

These take a rank-two contingency table. Marginal frequencies are formed from
row and column sums, and the independence reference is their outer product.
Use `op_plugin(freqs_dirichlet(y, prior_matrix))` for nonuniform cell priors.

## Tests

- `gstat(y, expected [, unit])`
- `chi2stat(y, expected [, unit])`
- `gstatindep(y [, unit])`
- `chi2statindep(y [, unit])`

Results are `type(test_result)` with real `stat`, integer `df`, and real `pval`.
Goodness-of-fit functions take vectors; independence tests take matrices.
Expected weights are normalized; sample counts determine sample size.

## Discretization

- `call discretize(x, num_bins, counts [, range, status])`
- `call discretize2d(x1, x2, num_bins1, num_bins2, counts [, range1, range2, status])`

`counts` is an allocatable integer vector/matrix output, allocated by the
procedure, including when already allocated by the caller. Ranges are real(dp)
length-two arrays. Bins are right-closed and include the lowest endpoint.
Optional status is zero on success or one on invalid inputs. Error outputs
are zero-filled with the requested shape; nonpositive bin counts produce a
zero-length extent. Input observations are never modified.

## Unsupported / intentionally different

`entropy.NSB` is untranslated because upstream contains only an external-process
interface, not its algorithm. No fake substitute or executable dependency exists.
Calling `entropy(..., method='NSB')` yields NaN. R attributes, labels, verbose
printing, recycled shapes, partial argument matching, and dynamic table/list
objects are not reproduced. Consult README.md for invalid-value, zero-cell,
precision and array-shape rules.
