# API coverage

Coverage basis: distinct exported computational functions in the upstream `NAMESPACE`. All 18 exported functions are computational and are included. No plotting, printing, data-download, or presentation-only exports are present.

Status: **substantial**. **18 of 18 (100%)** exported computational R functions have meaningful Fortran mappings.

The mapping count measures computational-function coverage, not exact R compatibility. R classes/attributes, RcppParallel execution, warning text, and exact online-worker floating-point ordering are outside the translated numerical API.

See the `Translation coverage` table in `README.md` and the machine-readable `[[extra.translation.function]]` entries in `fpm.toml`.
