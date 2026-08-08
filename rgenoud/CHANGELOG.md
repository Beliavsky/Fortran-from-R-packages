# Changelog

## 0.1.2

- Route the required lexical objective callback through a module-level typed dispatcher.
- Route the scalar objective used by the internal adapter through the existing typed dispatcher.
- Route the core cache evaluator through a typed module-level dispatcher as well.
- This removes the remaining direct callback invocations from nested internal procedures in the public adapter/cache paths, improving compatibility with gfortran versions that otherwise diagnose them as implicit interfaces under FPM.

## 0.1.1

- Fix gfortran/FPM portability for optional procedure callbacks in `rgenoud.f90`.
- Route optional analytical-gradient, lexical local-objective, local-gradient, and comparator callbacks through module-level procedures with non-optional typed procedure dummies.
- This avoids `-Werror=implicit-interface` failures in compilers that do not preserve the optional callback interface through nested internal-procedure host association.
- No numerical algorithm or public API change.

## 0.1.0

- Initial modern Fortran/FPM translation of rgenoud 5.9-0.3 computational code.
- Continuous and integer genetic optimization.
- Nine-operator accounting and operators 2-9.
- Scalar and lexical objective APIs.
- Self-contained BFGS/P9 local refinement.
- MemoryMatrix-style caching, numerical derivatives, Hessians, and statistics.
- Strict test suite and examples.
