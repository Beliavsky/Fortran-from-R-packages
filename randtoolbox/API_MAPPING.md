# API mapping

| R function / native kernel | Fortran counterpart | Notes |
|---|---|---|
| `setSeed` | `set_seed` | Convenience package seed; explicit seeds are preferred. |
| `congruRand` | `congru_rand`, `congru_rng` | General LCG, including 2^64 mode. |
| `SFMT` | `sfmt_generate`, `sfmt_rng` | All ten upstream Mersenne exponents and parameter sets. |
| `WELL` | `well_generate`, `well_rng`, `well_from_options` | All 17 rngWELL variants. |
| `knuthTAOCP` | `knuth_taocp`, `knuth_rng` | Knuth 2002 generator. |
| MT19937 used by `set.generator` | `mt19937_generate`, `mt19937_rng` | 2002 seed/array initialization and 32/53-bit resolution. |
| `torus` | `torus` | Pure and SFMT-mixed modes. |
| `halton`, `runif.halton` | `halton` | Pure and SFMT-mixed modes. |
| `sobol`, `runif.sobol`, `sobol.R` | `sobol` | Current upstream scrambling is disabled; dimensions <=1111. |
| `get.primes` | `get_primes`, `nth_prime` | Up to 100,000 primes. |
| `int2bit` | `int2bit` | Integer to bit vector. |
| `bit2int` | `bit2int` | Bit vector to integer. |
| `bit2unitreal` | `bit2unitreal` | Binary fraction conversion. |
| `%xor%` computational kernel | `bit_xor` | Bit-vector XOR. |
| `gap.test` | `gap_test` | Returns `rng_test_result`. |
| `freq.test` | `frequency_test` | Returns `rng_test_result`. |
| `serial.test` | `serial_test` | Returns `rng_test_result`. |
| `poker.test` / `doPokerTest` | `poker_test` | Native Fortran counting. |
| `order.test` | `order_test` | d=2..8. |
| `coll.test`, `coll.test.sparse` kernels | `collision_count`, `collision_test_counts` | Generator callback/UI layer separated from the collision computation. |
| `stirling` | `stirling_second` | Stirling numbers of the second kind. |
| `stirlingDividedByK` | `stirling_divided_by_k` | Translated recurrence. |
| `permut` | `permutations` | Recursive permutation construction. |
| `getWELLState` | `well_rng%get_state` | Explicit state object. |
| `rngWELLScriptR` | `well_rng` | Native recurrence replaces the R reference implementation. |
| `set.generator`, `put.description`, `get.description` | no direct R-global equivalent | R `.Random.seed`/callback integration; use explicit generator objects. |
| `trueRNG` helpers | omitted | External web-service access, not a numerical kernel. |
