# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of lmtest 0.9-40 computational code.
- Added OLS/WLS and generic coefficient/Wald/LR inference infrastructure.
- Translated Breusch-Godfrey, Breusch-Pagan, Durbin-Watson, Goldfeld-Quandt, Harvey-Collier, Harrison-McCabe, Rainbow, RESET, and Granger tests.
- Translated Cox, J, PE, and encompassing nonnested-model tests.
- Modernized the upstream Farebrother `pan.f` AS 153 routine as a free-format module.
- Added internal normal/Student-t/F/chi-square probability support.
- Added regression tests and example program.
- Retained upstream computational sources and metadata for provenance.
