# Computational coverage

## Translated

| Upstream area | Fortran implementation |
|---|---|
| `LL_HMM_Rcpp`, `nLL_hmm`, `ll_hmm` | `hmm_log_likelihood`, `forward_backward` |
| `LL_HHMM_Rcpp`, `nLL_hhmm` | `hhmm_log_likelihood`, `hhmm_forward_backward` |
| `viterbi`, `decode_states` | `viterbi_decode`, `decode_hhmm` |
| `simulate_hmm`, `simulate_observations` | `simulate_hmm_model`, `simulate_hhmm_model`, distribution RNGs |
| Distribution calculations for normal, lognormal, t, gamma, Poisson | `fhmm_distributions` |
| `Gamma2gammasUncon`, `gammasUncon2Gamma` and related transforms | `transition_to_unconstrained`, `unconstrained_to_transition` |
| `par2parUncon`, `parUncon2par` | HMM and HHMM `pack_*` / `unpack_*` routines |
| `get_initial_values` numerical heuristic | `initial_hmm_parameters`, `initial_hhmm_parameters` |
| `fit_model` numerical likelihood maximization | `fit_hmm`, `fit_hhmm` |
| `compute_ci` inverse-information calculations | covariance and standard-error fields in fit results |
| `compute_residuals` | `pseudo_residuals` |
| `predict.fHMM_model` | `forecast_hmm` |
| `reorder_states` | `reorder_hmm_states` |
| `compare_models`, AIC, BIC, logLik, npar, nobs | `compare_hmm_model` and fit-result fields |
| `compute_T_star` simulation branch | `compute_chunk_lengths` |

The Fortran implementation additionally exposes forward predictions, filtered
probabilities, and smoothed probabilities directly. These arrays are internal
to much HMM analysis but are not exposed as a compact standalone upstream API.

## Partially represented

- The R initializer uses repeated k-means and jittered candidate screening. The
  Fortran initializer uses sorted empirical quantiles, global scale estimates,
  persistent transition matrices, and optional multi-start jittering. The
  purpose is the same, but exact starting values differ.
- The upstream optimizer is R's `nlm`. The Fortran implementation uses native
  multi-start Nelder-Mead and numerical derivatives. Estimates need not be bit
  identical.
- `compute_T_star` date-calendar grouping is not compiled. Fixed and simulated
  chunk lengths are implemented; callers can supply empirically calculated
  `chunk_lengths` directly.
- Parameter fixing through strings such as `normal(mu=...)` is represented by
  callers constructing or packing only the parameterization they require; the
  R string parser is not translated.

## Excluded R infrastructure

- HTTP and market-data downloading.
- CSV/R data preparation built around data frames, dates, and `padr`.
- S3 print, summary, coefficient, and plotting methods.
- Colors, event annotations, and graphics.
- R parallel clusters and progress bars.
- RDS/RDA serialization and bundled fitted-model objects as compiled data.
- Formula, `xts`, and other R presentation adapters.

The unmodified source and bundled data remain under `original/` for provenance.
