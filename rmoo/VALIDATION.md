# Validation

The release was compiled with GNU Fortran 14.2.0 using:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

All permanent tests pass:

1. `test_pareto_reference` - nondominated fronts, crowding, Das-Dennis integer
   compositions/reference points, GD/IGD, nadir calculation, scalarizing and
   perpendicular-distance calculations.
2. `test_operators_survival` - SBX, polynomial mutation, HUX, OX, inversion,
   and NSGA-II/III/R-NSGA-II environmental selection.
3. `test_optimizers` - end-to-end real NSGA-II on ZDT1, real NSGA-III on DTLZ2,
   R-NSGA-II on ZDT2, plus binary, discrete-integer and permutation optimizers.
4. `test_upstream_regression` - reproduces the numerical portions of upstream
   rmoo tests: suggestions remain unchanged with zero crossover/mutation,
   NSGA-III result dimensions, perpendicular distance, and zero GD for an
   identical reference front.
5. `test_randomized_pareto` - 250 randomized 2-5 objective populations compared
   against an independent iterative Pareto-layer peeling implementation, with
   zero rank mismatches.

The ZDT1 example also builds and runs under the same flags.  In the release
validation it returned 45 nondominated members from a population of 60 after
50 generations.

No `.f` fixed-form files are present.  No `.f90` line exceeds 132 columns.
`fpm.toml` parses as TOML and explicitly sets free source form, no implicit
typing, and no implicit external interfaces.
