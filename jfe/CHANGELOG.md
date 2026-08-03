# Changelog

## 2.5.11-fortran.1

- Initial modern Fortran/FPM translation.
- Ported all 32 non-download computational exports.
- Added explicit status handling, IEEE-NaN omission, historical type-1 VaR,
  drawdown and CAPM utilities, four tests, three examples, and a demo.
- Preserved the upstream mixed decimal/percentage return conventions.
- Corrected the uninitialized modified Burke branch and honored `rf` in the
  nonannualized Sharpe ratio.
