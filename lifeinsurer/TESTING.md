# Testing

`run_tests.sh` builds with Fortran 2018 conformance, all common warnings as errors, bounds checking, backtraces, and floating-point traps. `run_release_tests.sh` repeats the test suite at `-O3 -Werror`.

The tests cover helper functions and rounding, all major cash-flow types, PVfactory recursions, premiums, prospective reserves, contract modification helpers, and profit participation.
