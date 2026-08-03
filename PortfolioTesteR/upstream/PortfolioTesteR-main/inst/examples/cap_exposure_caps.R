# ------------------------------------------------------------------------------
# Example: cap_exposure with per-symbol and per-group caps (and down-only renorm)
#
# Strategy outline:
#   1) Build a momentum strategy: each week, pick the top-10 by 12-month momentum
#      and equal-weight them (~10% each when 10 selected).
#   2) Apply exposure caps:
#        - Cap any single symbol at 9% (max_per_symbol = 0.09)
#        - Cap any “group/sector” at 45% total (max_per_group = 0.45)
#      Then *renormalize down only* so that if gross > 1 it’s scaled to 1.
#      Leftover (if any) is implicit cash (no leverage up).
#   3) Turn capped weights into portfolio returns and print basic metrics.
#
# Notes:
#   - We build a simple two-bucket group map (G1/G2) just for the demo.
#   - Uses weekly sample data bundled with the package.
#   - All calls use explicit `PortfolioTesteR::` to avoid masking surprises.
# ------------------------------------------------------------------------------

cat("\n=== cap_exposure: basic demo ===\n")                                      # banner so the example prints clearly

data(sample_prices_weekly, package = "PortfolioTesteR")                        # 1) load weekly adjusted prices shipped with the package
mom12 <- PortfolioTesteR::calc_momentum(sample_prices_weekly, 12)              # 2) compute 12-month momentum scores (panel, wide format)
sel10 <- PortfolioTesteR::filter_top_n(mom12, 10)                              # 3) select top-10 per date by those scores (wide logical selection)
w_eq  <- PortfolioTesteR::weight_equally(sel10)                                # 4) convert selection to equal weights (~0.1 each when 10 selected)
# Extract the universe
syms <- setdiff(names(w_eq), "Date")

# --- Build a demo group map (deterministic, reproducible, aligned to our universe) ---
# We need a stable "symbol → group" mapping that works the same on every machine.
# demo_sector_map() returns a two-column data.frame for exactly the symbols we pass
# (here: `syms`) and assigns them evenly into `n_groups` buckets "G1", "G2", ...
# This keeps the example reproducible while still showing per-group behaviour
# (e.g., sector-neutral caps or Top-K-per-group). To try more diversification,
# just change n_groups = 3L or 4L and re-run—no other code changes needed.
sector_map <- PortfolioTesteR::demo_sector_map(syms, n_groups = 2L)

# --- Convert the map to the most bullet-proof format for cap_exposure() ---
# cap_exposure() accepts either:
#   (A) a named character vector: names = symbols, values = group labels  ← safest across builds
#   (B) a data.frame with symbol/group columns (older builds were picky about lower-case)
# Using a *named vector* avoids column-name casing issues entirely and guarantees alignment:
#  - We coerce to character to avoid factor surprises on some R setups.
#  - We subset with [syms] so the vector matches the weight table’s columns 1:1 and
#    silently drops any extra rows that might be present in sector_map.
# If you really want to pass a data.frame, the new cap_exposure() normalises common
# column name variants, but the named-vector route keeps examples maximally portable.
group_map <- stats::setNames(as.character(sector_map$Group),
                             as.character(sector_map$Symbol))[syms]


w_cap <- PortfolioTesteR::cap_exposure(
  weights          = w_eq,
  max_per_symbol   = 0.09,
  group_map        = group_map,      # named vector → works on all builds
  max_per_group    = 0.45,
  renormalize_down = TRUE
)


pr <- PortfolioTesteR::portfolio_returns(w_cap, sample_prices_weekly)          # 8) turn weights into portfolio simple returns
m  <- PortfolioTesteR::perf_metrics(stats::na.omit(pr), freq = 52)             # 9) drop first NA from lag alignment; compute weekly metrics (52)

cat(sprintf("AnnReturn=%.2f%%  AnnVol=%.2f%%  Sharpe=%.2f  MaxDD=%.2f%%\n",    # 10) print concise summary
            100*m$ann_return, 100*m$ann_vol, m$sharpe, 100*m$max_drawdown))
