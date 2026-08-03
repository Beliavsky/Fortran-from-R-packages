# ------------------------------------------------------------------------------
# Example: limit_positions (legacy selection version) → equal-weight → metrics
# ------------------------------------------------------------------------------
#
# Strategy outline:
#   1) Build a 12-month momentum score panel.
#   2) Build a wide logical selection for the top-10 per date.
#   3) Enforce “no more than N positions per date” using the legacy
#      selection-based limit_positions() (keeps the highest scores when needed).
#   4) Convert the limited selection to equal weights (~1/N each) and print metrics.
#
# Why this example:
#   - Demonstrates the legacy (selection-based) limiter preserved for back-compat.
#   - Shows the English-like helpers `as_selection()` and `weight_equally()`.
# ------------------------------------------------------------------------------

cat("\n=== limit_positions (legacy selection) ===\n")                           # banner so the example prints clearly

data(sample_prices_weekly, package = "PortfolioTesteR")                        # 1) load weekly adjusted prices shipped with the package
mom12 <- PortfolioTesteR::calc_momentum(sample_prices_weekly, 12)              # 2) compute 12-month momentum scores (panel, wide format)

sel10    <- PortfolioTesteR::filter_top_n(mom12, 10)                           # 3) wide logical selection: top-10 per date
sel_mask <- PortfolioTesteR::as_selection(sel10)                               # 4) ensure selection is a clean wide logical table (Date + symbols)

lim <- PortfolioTesteR::limit_positions(                                       # 5) apply the legacy limiter
  selection_df   = sel_mask,                                                   #    (wide) logical selection
  max_positions  = 10L,                                                        #    hard cap per date (no more than 10)
  ranking_signal = mom12,                                                      #    tie-breaks/ranking by the momentum panel
  verbose        = FALSE                                                       #    keep output tidy
)

w_eq_lim <- PortfolioTesteR::weight_equally(lim)                               # 6) turn the limited selection into equal weights (~1/N each)

pr <- PortfolioTesteR::portfolio_returns(w_eq_lim, sample_prices_weekly)       # 7) portfolio simple returns from the weights
m  <- PortfolioTesteR::perf_metrics(stats::na.omit(pr), freq = 52)             # 8) drop first NA; compute weekly metrics (52)

cat(sprintf("AnnReturn=%.2f%%  AnnVol=%.2f%%  Sharpe=%.2f  MaxDD=%.2f%%\n",    # 9) print concise summary
            100*m$ann_return, 100*m$ann_vol, m$sharpe, 100*m$max_drawdown))

