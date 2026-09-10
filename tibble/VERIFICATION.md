# Verification

The translation has deterministic FPM tests for:

- heterogeneous construction and all four supported atomic column types;
- exact named extraction and explicit missing-value masks;
- validation of rectangular and mask-size invariants;
- row filtering and repeated/reordered row slicing;
- column selection, removal, insertion, replacement, and scalar recycling;
- schema-compatible row insertion and character-width preservation;
- homogeneous matrix conversion;
- typed `enframe()` and `deframe()` round trips;
- stable unique-name repair; and
- nonzero-row, zero-column tables.

Run:

```text
fpm test
```

During development the suite was also compiled with GNU Fortran runtime checking:

```text
fpm test --flag "-g -fcheck=all -fbacktrace"
```

The tests verify the documented Fortran API. They do not claim parity for the
R-runtime features listed as deliberately deferred in `API_COVERAGE.md`.
