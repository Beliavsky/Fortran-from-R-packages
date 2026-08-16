# API map

| R / C++ upstream | Fortran |
|---|---|
| `zrnorm(n)` | `zrnorm(n)` |
| `zrexp(n)` | `zrexp(n)` |
| `zrunif(n)` | `zrunif(n)` |
| `zsetseed(s)` | `zsetseed(s)` |
| `Ziggurat::rnorm()` | `rng%rnorm()` |
| `Ziggurat::rexp()` | `rng%rexp()` |
| `Ziggurat::runi()` | `rng%runi()` |
| `Ziggurat::kiss()` | `rng%kiss()` |
| `Ziggurat::setSeed()` | `rng%set_seed()` |
| `Ziggurat::getSeed()` | `rng%get_seed()` |
| `Ziggurat::getPars()` | `rng%get_state()` |
| `Ziggurat::setPars()` | `rng%set_state()` |

For array generation, the derived type also provides `fill_normal`,
`fill_exponential`, and `fill_uniform` to avoid temporary allocations.
