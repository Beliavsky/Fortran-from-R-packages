# stockAnalyst-fortran

A modern Fortran/FPM translation of the computational routines in the R package **stockAnalyst 1.0.1**.

The library implements all 55 exported equity-valuation formulas: dividend-discount models, free-cash-flow valuation, residual-income valuation, price and enterprise-value multiples, holding-period returns, growth estimates, CAPM/Fama-French required returns, and WACC.

## Build

```text
fpm build
fpm test
fpm run
fpm run --example reference_examples
```

No external numerical libraries are required. All public calculations use `real(dp)`, where `dp = kind(1.0d0)`.

## API style

Descriptive snake-case names are the preferred Fortran interface. Compatibility generic names matching the original R exports are also available; Fortran identifiers are case-insensitive. For example, both `computing_r_with_capm(...)` and `computingRwithCAPM(...)` resolve to the same implementation.

The original R package rounds every exported result. This port preserves those per-function decimal-place conventions using an explicit round-half-to-even routine.

See `API.md`, `PORTING.md`, and `TESTING.md` for details.

## License

GPL-3.0-only, matching the original package metadata. See `LICENSE` and `NOTICE.md`.
