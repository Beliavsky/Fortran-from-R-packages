# metaheuristicOpt-fortran

Modern Fortran translation of the computational code in the R package
`metaheuristicOpt` 2.0.0.

The original package implements 21 population-based continuous optimization
algorithms in pure R.  This project translates those numerical algorithms to
Fortran 2018 and provides an FPM project.  R progress bars, plotting, R list/data
frame handling, and documentation-generation machinery are intentionally not
translated.

## Algorithms

The public `metaopt()` dispatcher accepts the same algorithm names as the R
package:

| Name | Algorithm |
|---|---|
| `PSO` | Particle Swarm Optimization |
| `ALO` | Ant Lion Optimizer |
| `GWO` | Grey Wolf Optimizer |
| `DA` | Dragonfly Algorithm |
| `FFA` | Firefly Algorithm |
| `GA` | Genetic Algorithm |
| `GOA` | Grasshopper Optimisation Algorithm |
| `HS` / `IHS` | Harmony Search |
| `MFO` | Moth-Flame Optimizer |
| `SCA` | Sine-Cosine Algorithm |
| `WOA` | Whale Optimization Algorithm |
| `CLONALG` | Clonal Selection Algorithm |
| `DE` | Differential Evolution |
| `SFL` | Shuffled Frog Leaping |
| `CSO` | Cat Swarm Optimization |
| `ABC` | Artificial Bee Colony |
| `KH` | Krill-Herd Algorithm |
| `CS` | Cuckoo Search |
| `BA` | Bat Algorithm |
| `GBS` | Gravitational Based Search |
| `BHO` | Black Hole Optimization |

Each method is also available as its own Fortran subroutine (`pso`, `alo`, ...,
`bho`).  `metaopt_many()` runs several named algorithms with the same control
settings, corresponding to the vector-of-algorithm use of R's `metaOpt()`.

## Basic use

```fortran
program demo
   use metaheuristic_opt, only : dp, mh_control, mh_result, metaopt
   implicit none
   type(mh_control) :: control
   type(mh_result) :: result
   real(dp) :: lower(5), upper(5)

   lower = -10.0_dp
   upper =  10.0_dp
   control%num_population = 40
   control%max_iter = 200
   control%seed = 12345
   control%legacy_quirks = .false.

   call metaopt('GWO', sphere, lower, upper, result, control)
   print *, result%value
   print *, result%par
contains
   real(dp) function sphere(x) result(f)
      real(dp), intent(in) :: x(:)
      f = sum(x*x)
   end function sphere
end program demo
```

For maximization set `control%maximize = .true.`.

`mh_result` contains:

- `par`: best parameter vector returned by the selected method;
- `value`: objective value on the original minimization/maximization scale;
- `history`: per-generation best history where the source algorithm maintains one;
- `iterations`: number of generations executed;
- `evaluations`: objective callback count;
- `algorithm`: method name.

## Control defaults

`mh_control` starts with the defaults used by the package wrappers:

- common: `num_population=40`, `max_iter=500`;
- PSO: `vmax=2`, `ci=1.49445`, `cg=1.49445`, `w=0.729`;
- FFA: `b0=1`, `gamma=1`, `alpha_ffa=0.2`;
- GA: `pm=0.1`, `pc=0.8`;
- HS: `par=0.3`, `hmcr=0.95`, `bandwidth=0.05`;
- CLONALG: selection size `population/4`, multiplication factor `0.5`,
  hypermutation rate `0.1`;
- DE: scaling vector `0.8`, crossover rate `0.5`, strategy `"best 1"`;
- SFL: memeplex count `population/3`, frog-leaping iterations `10`;
- CSO: mixture ratio `0.5`, tracing constant `0.1`, maximum velocity `1`,
  `smp=20`, `srd=20`, `cdc=num_variables`, `spc=.true.`;
- ABC: cycle limit `num_variables*population`;
- KH: the eight numerical defaults in the original `KH()` wrapper;
- CS: abandoned fraction `0.5`;
- BA: maximum/minimum frequency `0.1/-0.1`, `gama=1`, `alpha_ba=0.1`;
- GBS: gravitational constant defaults to the maximum search bound and
  `kbest=0.1`.

DE accepts the source strategies `"clasical"`, `"best 1"`,
`"target to best"`, `"best 2"`, `"rand 2"`, and `"rand 2 dir"`.
For convenience the correctly spelled `"classical"` is accepted as an alias.

## Source-compatibility quirks

The R package contains a number of implementation details that are clearly
observable computational behavior.  They are available through
`control%legacy_quirks=.true.` (the default):

- PSO initializes `Gbest` through the package's inverted `calcBest()` call;
- FFA's final per-generation selector uses `which.max(Light)` although lower
  internal fitness is otherwise treated as better;
- CLONALG's final `calcBest()` call uses the opposite sign from most methods;
- SFL reshuffling sorts using the raw objective rather than `optimType`, which
  matters for maximization;
- DE and ABC preserve the source's cached fitness after post-selection bound
  clipping;
- CSO returns the best candidate recorded at the start of an iteration rather
  than forcing a final post-update rescan.

Set `legacy_quirks=.false.` for corrected versions of these behaviors.  The
algorithmic update formulas themselves are otherwise unchanged.

Some unsafe R-language artifacts are not reproduced literally.  For example,
Fortran performs explicit dimension/bound validation and does not rely on R's
partial recycling or non-integer indexing behavior.

## Package helper functions

The computational helpers from `metaheuristic.FunctionCollection.R` are also
provided:

- `sphere`
- `schwefel`
- `rastrigin`
- `cantilever_beam`
- `centilever_beam` (alias retaining the package's spelling)
- `mae`

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example basic_pso
```

The project has no external numerical-library dependency.

For validation of this release, FPM was not installed in the build environment,
so the same source layout was compiled directly with GNU Fortran 14.2.0 using
Fortran 2018 mode.

## Validation

The release test suite covers:

1. all 21 algorithms on a bounded four-dimensional sphere problem;
2. deterministic seeding and case-insensitive `metaopt()` dispatch;
3. maximization semantics;
4. all six DE strategies;
5. legacy/corrected FFA and CLONALG behavior;
6. odd-dimensional GOA;
7. the package benchmark/helper functions.

It has been run with both `-O2` and the bounds-checked configuration:

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Wimplicit-interface -Werror=implicit-interface
```

The library itself compiles without warnings under those flags.  Linkers may
print an executable-stack notice for test/example programs that pass internal
Fortran procedures as callbacks; that is a compiler trampoline detail, not a
library diagnostic.

## License and provenance

`metaheuristicOpt` 2.0.0 declares `GPL (>= 2) | file LICENSE`.  This translation
is therefore distributed under GPL-2.0-or-later.  The original package's
`LICENSE` is retained, the GPL v2 text is in `COPYING`, and the complete attached
R source tree is retained under `original/metaheuristicOpt-master/` for
provenance and algorithm auditing.
