# Validation

Compiler:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
```

Strict flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

`./scripts/test_gfortran.sh` passes all six regression programs:

1. `test_roi_regressions` - all three QPs from the upstream R package test file.
2. `test_qproblemb_hotstart` - simply bounded solve and persistent hotstart.
3. `test_general_constraints` - equality and two-sided general constraints.
4. `test_statuses` - infeasible QP and unbounded zero-Hessian LP.
5. `test_sqproblem_hotstart` - SQProblem-style warm matrix/gradient update.
6. `test_model_diagnostics` - objective/state/constraint-count getters.

Both examples also compile and run under the same flags.

Upstream ROI regression 1:

```text
x = 0.4761904762  1.0476190476  2.0952380952
objective = -2.3809523810
```

Upstream ROI regression 2:

```text
x = 2.0  0.0
objective = -2.0
```

Upstream ROI maximization regression:

```text
x = 0.4761904762  1.0476190476  2.0952380952
objective = 2.3809523810
```

The FPM manifest is syntactically valid TOML and all translated Fortran source
lines are at most 132 columns.  FPM itself is not installed in the validation
container, so the identical FPM source/test/example tree is compiled directly
with GNU Fortran using the strict scripts.

A final clean-archive extraction/rebuild is performed before release.
