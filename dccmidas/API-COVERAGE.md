# API coverage

The translation coverage basis is the set of distinct **exported computational R functions**. Plotting, printing, summary presentation, and other presentation-only methods are excluded. Internal helpers and operators that are not exported are also excluded.

Status: **substantial**

Coverage: **22 of 22 (100%)** mapped computational functions.

The 100% figure means that every function in the computational coverage basis has a meaningful Fortran implementation. It does **not** mean full compatibility with R lists, xts/zoo indexing, S3 objects, out-of-sample workflow code, optimizer trajectories, or presentation methods.

The authoritative per-function mapping is in `[extra.translation]` and `[[extra.translation.function]]` entries in `fpm.toml`. The same table is reproduced in `README.md`.
