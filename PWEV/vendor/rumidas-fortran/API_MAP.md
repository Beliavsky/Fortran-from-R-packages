# API coverage map

## Lag and data helpers

| R function | Fortran counterpart | Coverage |
|---|---|---|
| `beta_function` | `beta_function`, `beta_weights` | Direct |
| `exp_almon` | `exp_almon`, `exponential_almon_weights` | Direct |
| `mv_into_mat` | `mv_into_mat`, `lag_matrix_from_period_index` | Adapted: explicit integer period mapping replaces `xts` dates |
| `multi_step_ahead_pred` | `multi_step_ahead_pred` | Direct model recursion with explicit final state |

## GARCH-MIDAS models

Every function below is public with the same case-insensitive Fortran name.
The Fortran routines are subroutines returning allocatable arrays and an
explicit status value.

| R family | Fortran routines |
|---|---|
| GM skewed | `gm_loglik`, `gm_cond_vol`, `gm_long_run_vol` |
| GM non-skewed | `gm_loglik_no_skew`, `gm_cond_vol_no_skew`, `gm_long_run_vol_no_skew` |
| GM-X skewed | `gm_x_loglik`, `gm_x_cond_vol`, `gm_x_long_run_vol` |
| GM-X non-skewed | `gm_x_loglik_no_skew`, `gm_x_cond_vol_no_skew`, `gm_x_long_run_vol_no_skew` |
| GM-2M skewed | `gm_2m_loglik`, `gm_2m_cond_vol`, `gm_2m_long_run_vol` |
| GM-2M non-skewed | `gm_2m_loglik_no_skew`, `gm_2m_cond_vol_no_skew`, `gm_2m_long_run_vol_no_skew` |
| DAGM skewed | `dagm_loglik`, `dagm_cond_vol`, `dagm_long_run_vol` |
| DAGM non-skewed | `dagm_loglik_no_skew`, `dagm_cond_vol_no_skew`, `dagm_long_run_vol_no_skew` |
| DAGM-X skewed | `dagm_x_loglik`, `dagm_x_cond_vol`, `dagm_x_long_run_vol` |
| DAGM-X non-skewed | `dagm_x_loglik_no_skew`, `dagm_x_cond_vol_no_skew`, `dagm_x_long_run_vol_no_skew` |
| DAGM-2M skewed | `dagm_2m_loglik`, `dagm_2m_cond_vol`, `dagm_2m_long_run` |
| DAGM-2M non-skewed | `dagm_2m_loglik_no_skew`, `dagm_2m_cond_vol_no_skew`, `dagm_2m_long_run_no_skew` |

`garch_midas_evaluate` is the preferred typed interface and handles all these
families through `type(garch_midas_spec)`.

## MEM models

| R function family | Fortran counterpart | Coverage |
|---|---|---|
| `MEM_loglik`, `MEM_pred` | `mem_loglik`, `mem_pred` | Direct |
| non-skewed MEM | `mem_loglik_no_skew`, `mem_pred_no_skew` | Direct |
| MEM-X | `mem_x_loglik`, `mem_x_pred` and non-skewed forms | Direct |
| MEM-MIDAS | `mem_midas_loglik`, `mem_midas_pred`, `mem_midas_lr_pred` and non-skewed forms | Direct |
| MEM-MIDAS-X | `mem_midas_x_loglik`, `mem_midas_x_pred`, `mem_midas_x_lr_pred` and non-skewed forms | Direct |

`mem_evaluate` is the preferred typed interface and uses `type(mem_spec)`.

## Estimation and inference

| R function | Fortran counterpart | Coverage |
|---|---|---|
| `ugmfit` | `ugmfit`, `fit_garch_midas` | Adapted typed result; same model constraints and random-start idea |
| `umemfit` | `umemfit`, `fit_mem` | Adapted typed result |
| `Inf_criteria` | `information_criteria` | Direct |
| `LF_f` | `volatility_loss_functions` | Direct |
| `QMLE_sd` | robust covariance fields in `rumidas_fit_result` | Direct sandwich calculation from numerical observation scores |
| `MEM_QMLE_sd` | robust fields and MEM fallback scale | Adapted |

## Omitted R infrastructure

- `print.rumidas` and `summary.rumidas`;
- `xts`, `zoo`, `lubridate`, and calendar-frequency dispatch;
- packaged `.rda` data and documentation-only examples;
- automatic ARIMA estimation for an X-variable forecast;
- R list/data-frame formatting and captured console output.
