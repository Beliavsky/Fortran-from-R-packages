# adagio-fortran

Modern Fortran/FPM translation of the computational routines in the R package
**adagio 0.9.2** (Hans W. Borchers), a collection of discrete and global
optimization algorithms.

## Scope

The library translates the exported computational surface of `adagio`:

- derivative-free optimization: Nelder-Mead, bounded Nelder-Mead, and
  Hooke-Jeeves;
- global/stochastic optimization: `simpleEA`, `simpleDE`, and `pureCMAES`;
- numerical gradients;
- assignment, subset sum, 0/1 knapsack, multiple knapsack, set cover,
  change making, and approximate bin packing;
- Hamiltonian path/cycle search;
- maximum-sum subvector/submatrix and maximum empty rectangle;
- transfinite bound transforms;
- value counting and subsequence occurrence search;
- historized objective-value/input storage;
- the Rosenbrock, Rastrigin, Nesterov, Hald, Shor, Trefethen/Wagon and
  generated max-quadratic test functions.

The plotting-only `fminviz()` and `flineviz()` routines are intentionally not
translated.  R closures, `...`, list/S3 presentation and other language glue
are represented by explicit Fortran callbacks and derived types instead.

## lpSolve dependency

The original R package imports `lpSolve` for `assignment`, `change_making`,
`setcover`, and `mknapsack`.  This distribution vendors the user-supplied
**lpSolve-fortran v0.1.0** translation under
`vendor/lpSolve-fortran-v0.1.0/` and declares it as an FPM path dependency.
Those four routines formulate the same LP/MILP models as the R source; they do
not use unrelated replacement algorithms.

`adagio-fortran` is GPL-3.0-or-later.  The vendored lpSolve translation remains
LGPL-2.0-only; see `LICENSES.md` and the dependency's own `LICENSE`/`NOTICE.md`.

## Main API mapping

| R | Fortran |
|---|---|
| `neldermead` | `neldermead` |
| `neldermeadb` | `neldermeadb` |
| `hookejeeves` | `hookejeeves` |
| `simpleEA` | `simple_ea` |
| `simpleDE` | `simple_de` |
| `pureCMAES` | `pure_cmaes` |
| `assignment` | `assignment` |
| `subsetsum`, `sss_test` | `subsetsum`, `sss_test` |
| `knapsack`, `mknapsack` | `knapsack`, `mknapsack` |
| `setcover`, `change_making` | `setcover`, `change_making` |
| `bpp_approx` | `bpp_approx` |
| `maxsub`, `maxsub2d` | `maxsub`, `maxsub2d` |
| `maxempty` | `maxempty` |
| `hamiltonian` | `hamiltonian` |
| `transfinite` | `transfinite_forward`, `transfinite_inverse` |
| `count`, `occurs` | `count_values`, `occurs` |
| `Historize` | `history_buffer` |
| `maxquad` | `make_maxquad` + type-bound `value`/`gradient` |

All public symbols are re-exported from module `adagio`.

## Build with FPM

```text
fpm build
fpm test
fpm run --example optimization_demo
fpm run --example discrete_demo
```

The lpSolve dependency is local to the source tree, so no separate lpSolve
installation is required.

## Example

```fortran
program demo
  use adagio
  implicit none
  type(opt_result) :: fit

  fit = neldermead(fn_rosenbrock, [-1.2_dp, 1.0_dp])
  print *, fit%x, fit%f
end program demo
```

The included optimized example converges to approximately `(1,1)` with a
Rosenbrock objective around `7e-13` in the validation environment.

## Source-compatible quirks

The translation intentionally preserves observable numerical behavior from the
R source, including some unusual accounting/sampling details:

- `simpleDE` increments its reported `nfeval` by `N` for every trial point.
  `de_result%nfeval` preserves that value and `actual_nfeval` gives the real
  callback count.
- `simpleEA` retains the original population initialization and candidate-range
  indexing expressions, including their dependence on side lengths rather
  than an added lower-bound offset for random individuals.
- Hooke-Jeeves retains the package's target test based on `abs(fx) < target`.
- exact equality is retained where `count()` and `occurs()` use R's exact
  numeric comparison semantics.

See `TRANSLATION_NOTES.md` for details.

## Validation

Six regression executables cover the LP/MILP-backed discrete routines,
combinatorial routines, geometry/sequence utilities, test functions and
analytic gradients, numerical gradients/maxquad, local optimizers, stochastic
optimizers, CMA-ES, and historized storage.

The exact release sources pass with GNU Fortran 14.2.0 under both `-O2` and a
bounds-checked build using:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Werror
-Wimplicit-interface
```

GNU ld may print an executable-stack notice for executables that pass an
internal Fortran procedure as a callback; this is a compiler trampoline detail,
not a test failure.
