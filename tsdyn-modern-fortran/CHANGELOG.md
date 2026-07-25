# Changelog

## 0.1.0

- Added modern Fortran modules for AR, SETAR, LSTAR, local-linear AR, VAR, VECM, TVAR, and TVECM models.
- Added simulation, forecasting, model/order/rank selection, regime classification, IRFs, FEVD, GIRFs, bootstrap paths, rolling forecasts, and resampling.
- Added delta, BBC, Kapetanios-Shin, and Johansen numerical statistics.
- Added direct and box-index local-neighbor searches.
- Added demo, threshold example, and dated CSV fitting application.
- Added three runtime regression suites and strict debug/release build targets.
- Preserved GPL-2.0-or-later licensing in all Fortran files and added machine-enforced license checks.
- Fixed bounds defects found by multistep VAR forecast and bootstrap tests before release.
