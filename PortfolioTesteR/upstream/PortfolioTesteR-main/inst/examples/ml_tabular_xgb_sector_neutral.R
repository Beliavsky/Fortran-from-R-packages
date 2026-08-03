suppressPackageStartupMessages(library(PortfolioTesteR))
suppressPackageStartupMessages(library(data.table))

data(sample_prices_weekly); data(sample_prices_daily)

X <- ml_prepare_features(sample_prices_weekly, sample_prices_daily)
Y <- make_labels(sample_prices_weekly, horizon = 4L, type = "log")

symbols <- setdiff(names(sample_prices_weekly), "Date")
gmap <- demo_sector_map(symbols, n_groups = 4L)

xgb <- ml_make_model("xgboost",
                     params = list(objective="reg:squarederror", max_depth=3, eta=0.05),
                     nrounds = 200)

set.seed(2)
res <- ml_backtest(
  features_list = X, labels = Y,
  fit_fn = xgb$fit, predict_fn = xgb$predict,
  schedule = list(is = 104L, oos = 4L, step = 4L),
  group = "per_group", group_map = gmap,
  transform = "zscore",
  selection = list(top_k = 15L),
  weighting = list(method = "softmax", temperature = 12),
  prices = sample_prices_weekly,
  initial_capital = 1e5,
  name = "XGBoost Sector-Neutral"
)

print(res$backtest)
plot(res$backtest, type = "performance")
