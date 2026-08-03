# LifeInsureR-fortran

Modern Fortran 2018 implementation of the reusable computational core of the R package **LifeInsureR 1.0.1**.

The library models traditional life-insurance contracts through explicit typed inputs rather than R6 objects. It provides mortality and morbidity transitions, unit cash-flow generation, actuarial present values, premium and sum-insured calculations, prospective reserves, surrender and premium-free values, frequency corrections, costs, rounding, profit participation, contract extensions, premium waivers, and contract-grid premiums.

## Build

```sh
fpm test
fpm run
```

Direct GNU Fortran validation is available through `run_tests.sh`, `run_release_tests.sh`, and `run_tests.bat`.

## Scope

The compiled library covers the numerical contract engine. Excel export, RStudio project templates, R6/S3 object plumbing, HTML display, and automatic loading of `MortalityTables` objects are omitted. Mortality and morbidity probabilities are passed as ordinary arrays.

See `PORTING.md`, `API.md`, and `TRANSLATION_COVERAGE.md` for details.
