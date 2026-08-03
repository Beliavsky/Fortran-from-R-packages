# Build report

- Compiler: GNU Fortran 14.2.0
- Language mode: Fortran 2018
- Checked build: PASS
- Optimized build: PASS
- Tests: `test_analysis`, `test_portfolios`, and `test_prepare_data` all PASS
- Examples: `basic_portfolios` and `demo_ren` compile and run
- FPM manifests: parsed successfully as TOML for REN and vendored corpcor
- FPM executable: not installed in the validation environment; the package was
  therefore validated with the included source-order build scripts
