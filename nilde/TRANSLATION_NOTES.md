# Translation notes

## Scope

This release translates every exported computational routine in nilde 1.1-7:
`nlde`, `get.partitions`, `get.subsetsum`, `get.knapsack`, `bin.packing` and
`tsp_solver`.

S3 print methods, package-startup messages, R environments used as mutable
scratch storage, and `parallel` worker management are omitted.

## Diophantine engine

The R implementation writes the Voinov-Nikulin nested sums as recursive R
functions which enumerate the largest coefficient first and solve the smallest
coefficient algebraically. The Fortran implementation keeps the same search
structure but factors it into one reusable bounded Diophantine enumerator.
Coefficients are stably sorted internally and solutions are restored to the
original variable order.

The following modes are represented:

- nonnegative equality with at-most `M` parts;
- nonnegative equality with exactly `M` parts;
- 0/1 equality (`option=1`);
- 0/1 inequality (`option>1`) through the same slack-variable construction as
  the R package.

## Partitions

The R source uses a generating-function-derived recurrence. The Fortran routine
enumerates the mathematically identical canonical representation directly as a
nondecreasing vector of length `M`. For `at_most=.true.`, zero entries pad
partitions that contain fewer than `M` positive parts.

## Knapsack

The package reduces knapsack to a Diophantine equation by adding a unit slack
variable. The Fortran implementation enumerates the equivalent capacity
inequality directly, avoiding storage of the slack coordinate.

The original unbounded branch does not perform its advertised objective
filter. This is preserved by default and can be disabled with
`legacy_unbounded_all=.false.`.

## Bin packing

The original R code first enumerates every feasible 0/1 subset that fits in a
bin and then recursively combines these subsets. The Fortran translation uses
the equivalent exact set-partition search directly. Bin labels are canonical,
which removes pure bin-order permutations while retaining all distinct item
partitions.

## TSP

The package's exact TSP logic is preserved at the mathematical level:

1. compute an assignment relaxation lower bound;
2. compute a heuristic feasible upper bound;
3. increase the integer objective from the lower bound;
4. enumerate assignment-feasible selected arcs at that objective;
5. reject disconnected cycle covers/subtours;
6. stop at the first objective having one or more Hamiltonian tours.

The R dependencies used only in step 1/2 are replaced by native Hungarian and
cheapest-insertion routines. The exact search combines the original subset-sum
and degree-filter phases into a direct assignment-feasible recursion, which
avoids constructing enormous intermediate binary matrices but has the same
accepted tours.

## Integer precision

The R package stores integer values in R doubles. The Fortran port uses signed
64-bit integers for coefficients, right-hand sides, combinatorial solutions and
integer TSP costs. This is exact over the full int64 range subject to ordinary
arithmetic overflow limits.
