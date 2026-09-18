# stats test support

This development-only FPM package provides common assertions for the `stats`
translation tests. It is listed under `[dev-dependencies]`, so applications
using `stats-fortran` do not acquire test-support code as a production
dependency.

The assertion procedures accept absolute or relative tolerances and report the
maximum discrepancy before terminating a failed test. All current `stats`
tests use this module instead of defining local scalar, vector, matrix,
logical, or integer assertion procedures.

`migrate_test_assertions.py` performs the mechanical migration for newly
imported tests. It removes recognized local helpers, preserves host-defined
tolerances at call sites, wraps generated calls at 100 columns, and is
idempotent.
