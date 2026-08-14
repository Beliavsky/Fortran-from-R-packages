# API mapping

| rngWELL / C API | Fortran API |
| --- | --- |
| `setSeed4WELL` / MT2002 seed expansion | `well_rng%seed`, `init_mt2002` |
| `WELL2test` | `well_from_options` plus `%fill` / `%fill_matrix` |
| `InitWELLRNGXXX` | `well_rng%init('XXX', seed=...)` or `state=...` |
| `WELLRNGXXX` | `well_rng%next()` |
| raw WELL word | `well_rng%next_uint32()` |
| `GetWELLRNGXXX` | `well_rng%get_state()` |
| `putRngWELL` | `well_rng%put_state()` |
| `getRngWELL` | `%get_variant()` plus `%get_state()` |
| `initMT2002` | `init_mt2002()` |

The Fortran design uses explicit generator objects instead of global C state,
so multiple independent WELL streams can coexist safely in one process.
