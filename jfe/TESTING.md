# Testing

`run_tests.sh` builds with:

- Fortran 2018 conformance
- all common warnings and pedantic warnings
- warnings treated as errors
- array bounds and runtime checks
- floating-point traps for invalid operations, division by zero, and overflow

`run_release_tests.sh` builds with `-O3 -Werror` and runs every test, example,
and demo.

The four tests cover:

1. Annualization, VaR/ES, downside/upside moments, and Sharpe ratios.
2. Decimal and percentage drawdown paths, Calmar, Sterling, Pain, Ulcer, and
   modified Burke scaling.
3. Active premium, tracking error, information ratio, Jensen alpha, Treynor,
   appraisal, and M2 Sortino.
4. Tail and distributional ratios, adjusted Sharpe, annualized summary output,
   and Durbin h calibration.
