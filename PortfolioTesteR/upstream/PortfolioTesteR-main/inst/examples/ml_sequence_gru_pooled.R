suppressPackageStartupMessages(library(PortfolioTesteR))
suppressPackageStartupMessages(library(data.table))

data(sample_prices_weekly)

r1   <- panel_lag(calc_momentum(sample_prices_weekly,  1L), 1L)
r4   <- panel_lag(calc_momentum(sample_prices_weekly,  4L), 1L)
r12  <- panel_lag(calc_momentum(sample_prices_weekly, 12L), 1L)

features_seq <- list(r1 = r1, r4 = r4, r12 = r12)
Y <- make_labels(sample_prices_weekly, horizon = 4L, type = "log")

seq_gru <- ml_make_seq_model("gru", steps = 26L, units = 16L, epochs = 12L, seed = 42L)

set.seed(42)
res <- ml_backtest_seq(
  features_list   = features_seq,
  labels          = Y,
  steps           = 26L,
  horizon         = 4L,
  fit_fn          = seq_gru$fit,
  predict_fn      = seq_gru$predict,
  schedule        = list(is = 104L, oos = 4L, step = 4L),
  group           = "pooled",
  normalize       = "zscore",
  selection       = list(top_k = 10L),
  weighting       = list(method = "softmax", temperature = 12),
  prices          = sample_prices_weekly,
  initial_capital = 1e5,
  name            = "Seq-GRU: returns pyramid (26×3)"
)

print(res$backtest)
plot(res$backtest, type = "performance")
