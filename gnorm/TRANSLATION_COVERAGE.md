# Translation coverage

Upstream package: `gnorm` 1.0.2.

| R export | Fortran procedure | Status |
|---|---|---|
| `dgnorm` | `dgnorm` | Implemented |
| `pgnorm` | `pgnorm` | Implemented |
| `qgnorm` | `qgnorm` | Implemented |
| `rgnorm` | `rgnorm`, `rgnorm_fill` | Implemented |

Additional reusable procedures are provided for theoretical moments,
regularized gamma calculations, and explicit RNG state management.

There are no computational exports left untranslated. Plotting appears only in
examples and is not part of the Fortran library.
