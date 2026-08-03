suppressPackageStartupMessages(library(PortfolioTesteR))
suppressPackageStartupMessages(library(data.table))

data(sample_prices_weekly); data(sample_prices_daily)

# Build base features
X <- ml_prepare_features(sample_prices_weekly, sample_prices_daily)

# (Optional) sanity check to see available names
# print(names(X))  # should include: "mom12", "mom26", "mom52", "vol", "dist20", "rsi14"

# Add interactions using the actual names present in X
X <- ml_add_interactions(X, list(
  mom_vol = c("mom12","vol"),    # <— use "vol" here
  mom_rsi = c("mom12","rsi14")
))

Y <- make_labels(sample_prices_weekly, horizon = 4L, type = "log")

ridge <- ml_make_model("ridge")
rf    <- ml_make_model("rf", num.trees = 300, mtry = 3)
xgb   <- ml_make_model("xgboost",
                       params = list(objective = "reg:squarederror", max_depth = 4),
                       nrounds = 150)

ens <- ml_make_ensemble(ridge, rf, xgb)

set.seed(1)
res <- ml_backtest(
  features_list = X, labels = Y,
  fit_fn = ens$fit, predict_fn = ens$predict,
  schedule = list(is = 104L, oos = 4L, step = 4L),
  transform = "zscore",
  selection = list(top_k = 15L),
  weighting = list(method = "rank"),
  prices = sample_prices_weekly,
  initial_capital = 1e5,
  name = "Ensemble (Ridge+RF+XGB)"
)

print(res$backtest)
plot(res$backtest, type = "performance")

