# PortfolioTesteR 0.1.4

## New
- Add helpers for Chapter 3 workflow:
  - `ml_backtest_multi()` — run ML backtests across multiple horizons.
  - `pt_collect_results()` — gather coverage/IC/turnover/cost-sweep diagnostics.
  - `scores_oos_only()` — keep scores only on OOS decision dates.

## Docs
- Vignettes stabilized on knitr/rmarkdown (no Quarto dependency during check).

## Maintenance
- Extend globals to silence data.table NSE notes.
