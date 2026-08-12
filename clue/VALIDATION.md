# Validation

Compiler used for release validation:

```
GNU Fortran 14.2.0
```

## Debug/bounds-checked build

The library and tests were compiled with:

```
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

All regression executables passed.

The only link-time diagnostic on GNU/Linux is the executable-stack notice
caused by GNU Fortran trampolines for internal procedure callbacks used by the
SUMT convenience wrappers/tests.  It is not a compile-time or numerical test
failure.

## Optimized build

The complete suite was also compiled and run with:

```
-std=f2018 -O2
```

## Regression coverage

- `test_lsap_partition`: minimum/maximum LSAP, canonical labels, hard partition
  agreements, matching dissimilarities, lattice operations and fuzziness.
- `test_consensus_pclust`: DWH, hard/soft Euclidean and Manhattan consensus,
  including the LP-backed L1 stochastic update; fuzzy prototype clustering.
- `test_pava_medoid`: mean/median PAVA, exact LP/MILP k-medoids and PAM-style
  medoids.
- `test_transport_validity`: transportation-based Mallows distance and validity
  measures.
- `test_trees`: native ultrametric/additive-tree IP/IR kernels and
  ultrametrification.
- `test_target_fit`: fixed-topology L2/L1 ultrametric fits.
- `test_sumt`: generic sequential unconstrained minimization.
- `test_highlevel_fit`: high-level SUMT/L1/IRIP and sum-of-ultrametrics paths.

The final release archive is rebuilt and retested from a fresh extraction as a
packaging integrity check.


## v0.1.1 compiler-compatibility fix

`clue_sumt.f90` was restructured to avoid invoking optional host-associated procedure dummies from an internal procedure. The replacement uses module procedures with explicit abstract interfaces. All eight tests pass with `-Werror=implicit-interface` in both bounds-checked and optimized builds under GNU Fortran 14.2.0.
