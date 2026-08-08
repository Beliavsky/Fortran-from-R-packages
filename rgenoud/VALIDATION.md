# Validation

The translated sources were compiled with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Nine tests cover:

1. continuous Rosenbrock minimization with local refinement;
2. maximization;
3. exact integer optimization;
4. default lexical optimization;
5. custom lexical ordering;
6. analytical-gradient use and numerical Hessian output;
7. MemoryMatrix duplicate suppression;
8. boundary enforcement, P9 mixing, and crossover-count parity;
9. translated sample-moment diagnostics.

Both shipped examples were also compiled and run. The release archive is
validated a second time after extraction into a clean directory.

FPM was not installed in the validation container, so the same FPM source tree
was compiled directly with gfortran and `fpm.toml` was parsed independently.


## 0.1.1 portability regression

The public callback adapters in `src/rgenoud.f90` were additionally checked to
ensure that no optional procedure dummy is called directly from a nested internal
procedure. Optional callbacks are passed to module-level typed dispatchers, whose
procedure dummies are non-optional and have explicit abstract interfaces. This
addresses the Windows/gfortran FPM failure:

```text
Error: Procedure 'local_gradient' called with an implicit interface
```

The complete strict test suite was rerun after the change and again from a fresh
extraction of the release archive.

## 0.1.2 portability regression

The remaining required lexical objective callback (`fn`) and the core cached
evaluator (`evalfit`) are now routed through module-level typed dispatchers as
well. Consequently, no nested internal procedure in `src/rgenoud.f90` or the
cache adapter in `src/rgenoud_core.f90` directly invokes a host-associated user
procedure dummy. This addresses the subsequent Windows/gfortran diagnostic:

```text
Error: Procedure 'fn' called with an implicit interface
```

The full strict build and all nine tests were rerun after this broader fix.
