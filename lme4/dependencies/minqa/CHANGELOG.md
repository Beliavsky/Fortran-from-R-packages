# Changelog

## 0.1.0 - 2026-08-05

- Initial modern Fortran/FPM port of minqa 1.2.8.
- Added native APIs for BOBYQA, NEWUOA, and UOBYQA.
- Converted fixed-form numerical kernels to free-form module procedures.
- Added typed controls, results, workspace allocation, and status mapping.
- Added tests, example, documentation, and original-source provenance.
- Corrected the invalid `CALFUN(N,X,...)` call in the BOBYQA rescue path to
  evaluate the trial vector stored in `W`.
