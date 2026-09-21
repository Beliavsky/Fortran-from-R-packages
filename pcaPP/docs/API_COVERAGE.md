# API coverage

The coverage basis is **exported computational R functions**.  Plotting and
presentation exports are excluded.  In the supplied pcaPP NAMESPACE this gives
19 computational exports and three excluded presentation exports:
`PCdiagplot`, `objplot`, and `plotcov`.

All 19 computational exports have meaningful Fortran mappings.  The package
status is **substantial**, not complete, because the mapping fraction measures
translated computations rather than complete R-interface compatibility.

See the machine-readable `[[extra.translation.function]]` entries in
`fpm.toml` and the matching table in `README.md` for the exact mapping and
per-function status.
