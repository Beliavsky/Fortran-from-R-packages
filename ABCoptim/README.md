# ABCoptim-fortran

Modern Fortran translation of the computational core of the R package
**ABCoptim 0.15.0**, an implementation of the Artificial Bee Colony (ABC)
optimization method.

## Scope

The library translates both numerical implementations shipped by ABCoptim:

- `abc_optim`: behavior corresponding to the package's pure-R implementation,
  including its binary (`optiinteger`) mode.
- `abc_cpp`: behavior corresponding to the package's Rcpp implementation.

The translated algorithm contains the package's deterministic equally spaced
initial food sources, employed-bee phase, fitness-weighted onlooker phase,
greedy replacement, trial counters, one-scout-per-cycle rule, best-source
memory, and unchanged-best stopping criterion.

R S3 printing/plotting, `.Call`/Rcpp marshalling, graphics, and R-specific
object infrastructure are intentionally omitted.

## Build with FPM

```text
fpm build
fpm test
fpm run --example peaks
fpm run --example binary
```

The library has no external numerical dependencies.

## Basic API

```fortran
use abcoptim, only : dp, abc_control, abc_result, abc_cpp

type(abc_control) :: control
type(abc_result) :: result
real(dp) :: par(2), lb(1), ub(1)

par = 0.0_dp
lb = -10.0_dp
ub = 10.0_dp
control%seed = 213

call abc_cpp(par, objective, lb, ub, result, control)
```

Objective procedures have the interface

```fortran
function objective(x) result(value)
  real(dp), intent(in) :: x(:)
  real(dp) :: value
end function objective
```

Host association or module data can be used for additional objective data,
which replaces R's `...` mechanism.

## Controls

`abc_control` contains:

- `food_number = 20`
- `limit = 100`
- `max_cycle = 1000`
- `criter = 50`
- `optiinteger = .false.` (used by `abc_optim`)
- `seed`
- `legacy_r_scaling = .true.`

`parscale` and `fnscale` are optional arguments to both solver entry points.
As in ABCoptim, the optimized internal objective is
`fn(par/parscale)/fnscale`; a negative `fnscale` therefore turns minimization
into maximization.

## Result

`abc_result` includes the final food population, objective values, ABC fitness,
trial counters, best parameters/value, history, cycle count, actual objective
evaluation count, scout count, and a convergence flag. Fortran stores
`foods(dimension, food_source)` and `hist(dimension, history_index)` so each
point is contiguous in its first index.

## Compatibility details

The package's R and C++ implementations differ in several details. These are
kept separate rather than silently merged:

- The R implementation supports binary mutation; the C++ implementation does
  not.
- The R onlooker loop advances the source index after every trial and resets
  when it reaches `FoodNumber`, reproducing the package source behavior. The
  C++ loop advances only after a rejected probability test.
- R stops when persistence is greater than `criter`; C++ stops at greater than
  or equal to `criter`.
- Their probability formulas use slightly different tiny constants.
- C++ history includes the initial best source; R history begins with the first
  completed ABC cycle.
- The R implementation initializes `GlobalMin` with unscaled `fn(par)` even
  though population values use the scaled objective. This historical behavior
  is the default. Set `control%legacy_r_scaling=.false.` to initialize it with
  the scaled objective consistently.

Two unsafe wrapper/source edge cases are intentionally corrected: bounds may
be scalar or exactly dimension-sized without R's accidental vector
replication, and exhausting `max_cycle` returns a valid trimmed history rather
than reproducing the upstream off-by-one indexing at the hard limit.

The RNG is a deterministic native Fortran Park-Miller generator, not R's RNG.
Identical integer seeds therefore do not produce the same random trajectory as
`set.seed()` in R.

## Validation

The regression suite covers:

- ABC fitness transformation.
- The documented two-dimensional minimum at `(pi, pi)`.
- The package test function with minimizer near `-15.81515`.
- R-style binary optimization.
- `fnscale < 0` maximization.
- Deterministic equally spaced initial food sources.
- The R/C++ history convention difference.

See `TRANSLATION_NOTES.md` for a detailed source mapping.

## License

MIT. See `LICENSE` and `LICENSES.md`. The original package is retained under
`original/ABCoptim-master/` for provenance.
