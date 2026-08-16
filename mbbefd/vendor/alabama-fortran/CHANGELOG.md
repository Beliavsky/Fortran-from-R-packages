# Changelog

## 0.1.1

- Fixed GNU Fortran/FPM portability for optional procedure callbacks.
- Objective, gradient, constraint, and constraint-Jacobian invocations now use
  module-level explicit-interface trampolines instead of direct calls from
  internal procedures.
- Verified compilation with `-Werror=implicit-interface`.
- Removed generated object/module/build files from the source release.

## 0.1.0

- Initial modern Fortran computational translation of alabama 2025.1.0.
- Added `auglag`, `auglag1`, `auglag2`, `auglag3`, `constr_optim_nl`,
  `adpbar`, `augpen`, and `alabama_legacy`.
- Bundled supplied `numDeriv-fortran` as a local dependency.
- Bundled `roptim` as the native `optim`-style inner solver.
- Added typed controls/results, numerical KKT checks, tests and examples.
