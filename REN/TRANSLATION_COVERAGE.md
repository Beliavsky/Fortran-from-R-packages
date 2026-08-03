# Translation coverage

## Included

- All exported numerical portfolio constructors in `R/po_fun_test.R`.
- Hierarchical canonical-correlation clustering.
- Date filtering, month indexing, and missing-column removal.
- Rolling portfolio weights, turnover, daily gross returns, cumulative returns,
  wealth indices, annualized Sharpe ratios and volatility, and maximum drawdown.
- The intended covariance-shrinkage/Ledoit-Wolf slot in the rolling analysis.

## Excluded

- ggplot2 plots and reshape2 data preparation used only by plots.
- `setup_parallel`, interactive input, R cluster objects, foreach backends, and
  progress messages.
- Bundled `FF25.rda`; users pass typed arrays directly.

## Reused supplied translations

- corpcor for shrinkage, pseudoinverses, SVD, and positive-definite repair.

The supplied GPL-2.0-only glmnet translation is not linked because it is not
license-compatible with an AGPL-3.0-or-later combined work. The required
Gaussian elastic-net subset is implemented locally.
