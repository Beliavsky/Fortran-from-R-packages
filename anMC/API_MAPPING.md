# API mapping

Upstream package: `anMC` 0.2.5.

| R export | Fortran counterpart | Notes |
|---|---|---|
| `ANMC_Gauss` | `anmc_gauss` | Returns `mc_result`; full diagnostic data are available without a `typeReturn` switch. |
| `MC_Gauss` | `mc_gauss` | Supports optional pre-calibrated `mc_params`, matching upstream `params`. |
| `ProbaMax` | `proba_max` | Fixed `q` or `q_limits`; `ANMC`/`GANMC` and `MC`/`GMC` algorithm names. |
| `ProbaMin` | `proba_min` | Same structure as `proba_max`. |
| `conservativeEstimate` | `conservative_estimate` | Returns `conservative_result` with logical set membership, level, probability, and upstream-style uncertainty field. |
| `get_chronotime` | `chronotime_ns` | Returns a `real(dp)` nanosecond counter based on `system_clock`. |
| `mvrnormArma` | `mvrnorm_arma` | Result orientation is dimension x number-of-draws, like the C++ routine. `chol/=0` means an upper Cholesky factor is supplied. |
| `selectActiveDims` | `select_active_dims` | Implements methods 0, 1, 2, 3, 4, and 5. |
| `selectQdims` | `select_q_dims` | Implements the sequential q-increase stopping rule and the default q cap of 300. |
| `trmvrnorm_rej_cpp` | `trmvrnorm_rej_cpp` | Rejection sampler with adaptive batch size and explicit safety limits. |

## R concepts replaced by Fortran types

- R lists returned by `ProbaMax`/`ProbaMin` -> `probability_estimate`.
- R lists returned by `ANMC_Gauss`/`MC_Gauss` -> `mc_result`.
- R `problem` list -> `anmc_problem`.
- `attr(p, "error")` from `pmvnorm` -> `probability_result%error` internally and `pq_error` in the top-level result.
- R `lightReturn`/`typeReturn` switches -> one typed result containing the heavy information when it is computed.
- `set.seed()` -> `seed_fortran_rng`.

## Interface differences

The R package permits arbitrary R functions in `pmvnorm_usr` and `trmvrnorm`.
The Fortran v0.1.0 API exposes its native probability and truncated-normal
routines directly rather than carrying dynamic R callback semantics into the
core library.  The lower-level modules are public and can be adapted to a
project-specific procedure interface if a different backend is required.
