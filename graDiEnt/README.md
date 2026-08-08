# graDiEnt-fortran

Modern Fortran translation of the computational code in the R package
`graDiEnt` 1.0.1 by Brendan Matthew Galdo.

The package implements Stochastic Quasi-Gradient Differential Evolution
(SQG-DE), including the `rand/1/bin`, `current/1/bin`, and `best/1/bin`
adaptation schemes used by the R package.

## Build

```text
fpm build
fpm test
```

The library is standalone and has no external numerical dependencies.

## Main API

```fortran
use gradient

type(sqgde_options) :: opt
type(sqgde_result) :: res

call get_algo_params(5, opt)
opt%n_iter = 500
opt%adapt_scheme = SQGDE_RAND
call seed_rng(12345)
call optim_sqgde(my_objective, opt, res)
```

An objective has the explicit interface

```fortran
function my_objective(x) result(f)
  use gradient_kinds, only : dp
  real(dp), intent(in) :: x(:)
  real(dp) :: f
end function
```

`opt%init_center` and `opt%init_sd` are length-`n_params` arrays.  Fortran
scalar assignment can be used to broadcast a common value to all elements.

## Implemented computational features

- population initialization from independent normal distributions;
- SQG-DE stochastic quasi-gradient construction from `n_diff` vector pairs;
- `rand/1/bin`, `current/1/bin`, and `best/1/bin` base-vector schemes;
- binomial coordinate crossover with at least one selected coordinate;
- self-scaling `psi` normalization;
- additive uniform jitter;
- greedy minimization acceptance;
- optional objective re-evaluation (`purify`) for noisy objectives;
- standard-deviation and percent-improvement stopping criteria;
- thinned particle and objective traces;
- deterministic seeding through `seed_rng`.

R cluster orchestration (`FORK`, `PSOCK`, `doParallel`) is not translated.
The R implementation's parallel path is algorithmically synchronous, and the
Fortran serial implementation preserves that synchronous population update.

See `TRANSLATION_COVERAGE.md` for exact source-fidelity notes and defensive
edge-case fixes.

## License and provenance

The upstream package is distributed under the MIT license.  The original
package source is preserved under `original/graDiEnt-master/`.
