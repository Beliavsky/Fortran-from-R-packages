# psoptim-fortran 0.1.0

Modern Fortran translation of the computational code in the R package
`psoptim` 1.0 by Krzysztof Ciupke.

`psoptim` implements a basic particle swarm optimizer for **maximization**.
The original package is pure R; plotting/animation is deliberately omitted.

## Features

- Basic global-best particle swarm optimization.
- Inertia (`w`), cognitive (`c1`), and social (`c2`) terms.
- Per-coordinate position and velocity limits.
- Personal-best and global-best tracking.
- Mean- and best-fitness histories.
- User-selectable deterministic RNG seed.
- Final population, velocity, and personal-best populations returned for
  diagnostics.
- Literal compatibility switches for three quirks in `psoptim` 1.0.

## Compatibility details

The original R source contains three behaviors that are easy to overlook:

1. Initial positions in *every* dimension are sampled from
   `xmin[1] .. xmax[1]`, not each dimension's own bounds.
2. `runif(n*d, min=-vmax, max=vmax)` recycles the `vmax` vector over the
   linear random draw stream before the result is reshaped into a matrix.
3. If velocity is either above `+vmax[j]` **or below `-vmax[j]`, it is set
   to `+vmax[j]`**.

These are preserved by default with:

```fortran
control%legacy_initial_bounds = .true.
control%legacy_velocity_initialization = .true.
control%legacy_velocity_clip = .true.
```

For conventional per-dimension initialization and symmetric velocity
clipping, set all three to `.false.`.

The original R code repeatedly reevaluates `FUN(x.best.czastki)` and the
single current global-best point. The Fortran API accepts a pointwise
objective and caches those already known fitness values. This does not alter
PSO trajectories for ordinary deterministic pointwise fitness functions, but
it intentionally does not emulate side effects or population-dependent R
fitness functions.

## API

```fortran
use psoptim

type(ps_control) :: control
type(ps_result)  :: result

call ps_optimize(objective, xmin, xmax, vmax, result, control)
```

The callback is:

```fortran
function objective(x) result(f)
   use psoptim, only : dp
   real(dp), intent(in) :: x(:)
   real(dp) :: f
end function objective
```

The optimizer maximizes `f`, matching `psoptim`.

## Build with FPM

```text
fpm build
fpm test
fpm run --example rastrigin
```

The library has no external dependencies.

## Validation

All four regression executables pass under both optimized and bounds-checked builds.
The regression suite covers:

- two-dimensional Rastrigin maximization;
- the original dimension-1 initialization behavior and corrected mode;
- corrected multidimensional bounds and velocity clipping;
- best/mean histories and objective-evaluation accounting;
- literal negative-velocity clipping and recycled velocity initialization.

## Licensing

The original package declares `GPL (>= 2.0)`. This translation is therefore
released under GPL-2.0-or-later. See `COPYING`, `LICENSES.md`, and
`original/`.
