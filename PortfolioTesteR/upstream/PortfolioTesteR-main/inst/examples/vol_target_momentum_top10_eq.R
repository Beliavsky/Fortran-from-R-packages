# Vol-targeted momentum (top-10, equal-weight) — weekly sample
# Requires the patched PortfolioTesteR with:
# - vol_target() (cap=TRUE adds CASH)
# - portfolio_returns() (handles CASH even if prices don’t have it)

suppressPackageStartupMessages(library(PortfolioTesteR))

# --- Parameters you can tweak ---
lookback         <- 26     # weeks for cov/vol estimate
target_annual    <- 0.12   # target annualized vol (12%)
periods_per_year <- 52
cap_down_only    <- TRUE   # TRUE: scale down only (adds CASH)

# --- Data ---
data(sample_prices_weekly, package = "PortfolioTesteR")

# --- Simple momentum -> select top 10 -> equal-weight ---
mom12 <- PortfolioTesteR::calc_momentum(sample_prices_weekly, 12)
sel10 <- PortfolioTesteR::filter_top_n(mom12, 10)
w_eq  <- PortfolioTesteR::weight_equally(sel10)

# --- Vol targeting (adds CASH when cap_down_only = TRUE) ---
w_vt <- PortfolioTesteR::vol_target(
  weights          = w_eq,
  prices           = sample_prices_weekly,
  lookback         = lookback,
  target_annual    = target_annual,
  periods_per_year = periods_per_year,
  cap              = cap_down_only
)

# --- Portfolio returns (weekly) ---
pr_base <- PortfolioTesteR::portfolio_returns(w_eq,  sample_prices_weekly)
pr_vt   <- PortfolioTesteR::portfolio_returns(w_vt, sample_prices_weekly)

# drop any leading NA in returns
pr_base <- pr_base[is.finite(pr_base$ret), ]
pr_vt   <- pr_vt[is.finite(pr_vt$ret), ]

# --- Performance summary (weekly freq → 52) ---
pm_base <- PortfolioTesteR::perf_metrics(pr_base, freq = 52)
pm_vt   <- PortfolioTesteR::perf_metrics(pr_vt,   freq = 52)

cat("\n=== Momentum Top-10 EQ (weekly) ===\n")
cat(sprintf("Baseline: ann_return=%.2f%% ann_vol=%.2f%% Sharpe=%.2f max_dd=%.2f%%\n",
            100*pm_base$ann_return, 100*pm_base$ann_vol, pm_base$sharpe, 100*pm_base$max_drawdown))
cat(sprintf("VolTarget: ann_return=%.2f%% ann_vol=%.2f%% Sharpe=%.2f max_dd=%.2f%%\n",
            100*pm_vt$ann_return, 100*pm_vt$ann_vol, pm_vt$sharpe, 100*pm_vt$max_drawdown))

# --- Quick equity curve plot (base-R) ---
ec_base <- cumprod(1 + pr_base$ret)
ec_vt   <- cumprod(1 + pr_vt$ret)

op <- par(no.readonly = TRUE); on.exit(par(op))
par(mar = c(4, 4, 2, 1))
plot(pr_base$Date, ec_base, type = "l", xlab = "Date", ylab = "Equity (rebased to 1)",
     main = "Momentum Top-10 EQ vs Vol-Target (weekly)")
lines(pr_vt$Date,   ec_vt)
legend("topleft", bty = "n",
       legend = c("Baseline", "Vol-target"),
       lty = c(1,1))
