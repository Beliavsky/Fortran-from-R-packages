# Origin and provenance

- Original package: `tseries`
- Original version: 0.10-62
- Package title: Time Series Analysis and Computational Finance
- Original authors: Adrian Trapletti and Kurt Hornik
- BDS contributor: Blake LeBaron
- Source supplied as: `tseries-master.zip`
- Source archive SHA-256: `0d9b76941bad9ac0cff9abd335d8a311d2c96b6c984dc9316f38b7077f4756b9`
- Original package license declaration: `GPL-2 | GPL-3`
- Fortran translation date: 2026-07-22
- Translation status: experimental and independently produced

The original source contained R, C, and legacy Fortran code. This project reimplements the computational behavior in free-form modern Fortran and does not require the R runtime.

Omitted intentionally:

- Plotting and graphics methods
- Internet data download
- R S3 object presentation methods
- Most irregular-time-series, date, and file-format plumbing
- Bundled R datasets

The Fortran source contains attribution comments at module level through this project-wide provenance record. Original copyright notices remain credited in `NOTICE`.
