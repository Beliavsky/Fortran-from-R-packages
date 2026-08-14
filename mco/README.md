# mco-fortran

A modern Fortran translation of the computational core of the R package
`mco` 1.17.  The project is an FPM package and is licensed under
GPL-2.0-only, matching the upstream package.

## Implemented

- NSGA-II for bounded multiple-objective minimization
- Nonlinear constraints, feasible when every returned value is nonnegative
- Deb constrained dominance, nondominated sorting, crowding distance,
  binary tournaments, simulated-binary crossover, polynomial mutation,
  and elitist parent/offspring replacement
- Pareto masks and filtered fronts
- Front normalization, generational distance, generalized spread,
  dominated hypervolume, and the upstream additive epsilon indicator
- The exported Belegundu, Binh, Deb, Fonseca, Gianna, Hanne, Jimenez,
  VNT, and ZDT benchmark functions

Fortran matrices use `variables/objectives by points`: `par(nvar,npoint)`
and `value(nobj,npoint)`. This is the transpose of the row-oriented matrices
returned by the R interface.

## Build

```text
fpm build
fpm test
fpm run --example example_binh1
fpm run demo_mco
```

GNU Fortran scripts are also provided:

```text
./scripts/build_all.sh checked
./scripts/build_all.sh optimized
```

On Windows, run `scripts\\build_all.bat checked` or `optimized` from a
GNU Fortran command prompt.

## Minimal use

```fortran
use mco, only : dp, nsga2_options, nsga2_result, nsga2_optimize

type(nsga2_options) :: options
type(nsga2_result) :: result

options%population_size = 100
options%generations = 200
call nsga2_optimize(objective, 2, 2, [0.0_dp, 0.0_dp], &
                    [1.0_dp, 1.0_dp], result, options)
```

Objective callbacks fill `f(:)` and all objectives are minimized. Constraint
callbacks fill `g(:)` and are feasible when `g >= 0` elementwise.

See `doc/API_MAP.md` and `doc/PORTING_NOTES.md` for differences from R.
