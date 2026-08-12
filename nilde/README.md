# nilde-fortran

Modern Fortran 2018 translation of the computational core of the R package
**nilde 1.1-7**, *Nonnegative Integer Solutions of Linear Diophantine Equations
with Applications*.

The original package is by Natalya Pya Arnqvist, Vassilly Voinov, Rashid
Makarov and Yevgeniy Voinov and is distributed under GPL-2.0-or-later.
The original source tree is retained under `original/nilde-master/`.

## Implemented computational API

```fortran
use nilde
```

The public routines are:

- `nlde` -- enumerate nonnegative integer solutions of a linear Diophantine
  equation, including the package's binary equality and binary-inequality
  modes.
- `get_partitions` -- additive partitions into at most or exactly `M` parts.
- `get_subsetsum` -- 0/1 and bounded subset-sum enumeration.
- `get_knapsack` -- unbounded, 0/1 and bounded knapsack.
- `bin_packing` -- exact one-dimensional bin packing with all canonical optimal
  packings.
- `tsp_solver` -- exact integer-cost TSP using an assignment lower bound,
  heuristic upper bound, increasing objective levels, degree constraints and
  exact subtour rejection.
- `assignment_lower_bound` -- minimum assignment cost used by `tsp_solver`.

All integer combinatorial data use `integer(int64)` through the exported `i8`
kind. Objective coefficients use `real(real64)` through `dp`.

## Example

```fortran
use nilde

type(integer_solutions_t) :: r

r = nlde([3_i8,2_i8,5_i8,16_i8], 18_i8, m=6, at_most=.false.)
print *, r%nsol
print *, r%x
```

For the package documentation example this returns three exactly-six-part
solutions.

## Package-specific compatibility behavior

`nilde::get.knapsack(..., problem="uknap")` contains an observable historical
quirk: the R code adds a slack variable and returns **all feasible unbounded
vectors**, without filtering them to the maximum objective value. This is
preserved by default:

```fortran
r = get_knapsack(obj, a, cap, problem='uknap')
```

Set

```fortran
r = get_knapsack(obj, a, cap, problem='uknap', &
                 legacy_unbounded_all=.false.)
```

to obtain only objective maximizers. `r%objective` reports the best value in
both modes.

## TSP translation

The R package uses two optional R packages only for initialization:
`lpSolve::lp.assign()` supplies the lower bound and `TSP::solve_TSP()` supplies
a heuristic upper bound. The Fortran translation is self-contained:

- an O(n^3) Hungarian assignment routine supplies the lower bound;
- deterministic cheapest insertion supplies the default upper bound;
- `method='nearest_neighbor'` selects the alternate built-in heuristic;
- the exact search then follows the same objective-level logic as nilde:
  assignment-feasible edge selections are generated for each integer cost,
  degree constraints are enforced, and only a single Hamiltonian cycle is
  accepted.

The input cost matrix is integer-valued. Diagonal entries are ignored. A
negative off-diagonal entry means that the edge is unavailable; this is the
Fortran equivalent of an R `NA` edge.

The R `cluster` argument is not reproduced: parallel R worker orchestration is
interface infrastructure, not part of the numerical algorithm.

## Bin-packing result representation

The R implementation serializes bins into strings. The Fortran result instead
uses the computationally useful canonical matrix
`result%assignment(item,solution)`, where bin labels start at 1 in order of
first creation. `bin_ineff` and `total_ineff` reproduce the package's
inefficiency quantities.

## Build

```text
fpm build
fpm test
```

For an aggressive GNU Fortran debug build:

```text
fpm test --flag "-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Werror -Wimplicit-interface -Werror=implicit-interface"
```

FPM was not available in the translation environment, so the FPM source layout
was validated with equivalent direct GNU Fortran 14.2.0 builds.
