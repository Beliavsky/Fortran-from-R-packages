# Translation coverage

## Directly translated package behavior

The following `nls2` functionality is represented directly at the numerical-array level:

- brute-force/grid starting-value evaluation;
- random starting-value evaluation;
- Latin-hypercube starting-value evaluation;
- multiple explicit starting rows;
- choosing the minimum-RSS candidate;
- retaining all candidates/results;
- default multi-start nonlinear least squares;
- partially-linear brute-force/random/LHS evaluation;
- partially-linear nonlinear least squares;
- weights;
- singular-Jacobian start-point evaluation;
- the source grid-size rule and generated-start semantics.

## Standalone replacement for `stats::nls`

`nls2` delegates actual nonlinear optimization to R's `stats::nls`. To make the FPM package standalone, this translation supplies a Gauss-Newton solver with the same broad model/increment/step-halving structure and an `nls.control`-like convergence interface.

This is not a line-for-line port of every R `nls` algorithm:

- `default`: implemented as Gauss-Newton with step halving and QR-style convergence.
- `port`: implemented as the same core with bound projection; it is **not** a port of R's PORT/NL2SOL implementation.
- `plinear`: implemented with variable projection, solving the linear coefficients exactly and numerically differentiating the projected model with respect to nonlinear coefficients.

## Formula/model-object infrastructure omitted

The Fortran API deliberately does not reproduce:

- R formula parsing/evaluation;
- model frames, environments, `subset`, or `na.action` objects;
- S3 `nls`, `summary.nls`, `predict.nls`, `anova.nls`, printing, or updating methods;
- `proto` infrastructure.

Models are typed callbacks over numeric arrays.

## External suggested packages

### `lhs`

The R package calls `lhs::randomLHS`. The Fortran release implements ordinary randomized Latin-hypercube sampling directly and therefore has no external LHS dependency. The distributional construction is equivalent, but the random stream is not bit-identical to R/lhs.

### `CPoptim`

`nls.CPoptim` is only a wrapper around the separate `CPoptim` package; the convex-partition algorithm itself is not contained in the uploaded `nls2` source. The Fortran procedure `cpoptim_compat` provides bounded Latin-hypercube evaluation as a standalone compatibility fallback. It is **not** represented as an exact CPoptim translation.

## Source/documentation discrepancies

The current `nls2.R` source has two noteworthy discrepancies:

1. For explicit `start` tables with more than two rows, the help page says `random-search` samples rows without replacement, but the current source evaluates the supplied rows verbatim. The Fortran translation follows the current source.
2. The current source tests for `"plinear-brute"` in one two-row branch even though the accepted full algorithm name is `"plinear-brute-force"`. The Fortran translation implements the documented/intended grid behavior for both spellings instead of preserving this apparent typo.

## Numerical differences

- Fortran QR uses reorthogonalized modified Gram-Schmidt rather than R's LAPACK-backed QR.
- Numerical derivatives and RNG streams will not be bit-identical to R.
- Exact zero-residual problems may require nonzero `scale_offset`, as with R `nls`.
