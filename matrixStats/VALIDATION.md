# Validation record

Validated on 2026-10-11 with GNU Fortran 14.2.0, FPM 0.13.0, R 4.5.0 and the
pinned matrixStats 1.5.0 source described in NOTICE.md. Execution platform was
Linux. The manifest and sources avoid platform-specific system-library links;
native Windows execution was not available and is not claimed.

## Final build/test/example checks

An independent cleanly scoped check ran all nine stages against unchanged
package contents. Every stage passed:

| Profile | Build | Tests | All examples |
|---|---|---|---|
| Default | Pass | Pass | Pass |
| Debug with `-fcheck=all -fbacktrace` | Pass | Pass | Pass |
| Release | Pass | Pass | Pass |

Reproduce from the package root:

```text
fpm build
fpm test
fpm run --example
fpm build --profile debug --flag "-fcheck=all -fbacktrace"
fpm test --profile debug --flag "-fcheck=all -fbacktrace"
fpm run --example --profile debug --flag "-fcheck=all -fbacktrace"
fpm build --profile release
fpm test --profile release
fpm run --example --profile release
```

No fast-math or finite-only optimization was used. The final suite has **2,694
checks** in each profile:

- 237 deterministic end-to-end R-reference comparisons, covering all 99 mapped
  computational exports, selected/missing variants and output orientation
- 2,335 deterministic R ordering/weighted oracle checks: quantile types 1–9,
  six deterministic rank tie methods, all weighted median policies, unequal
  duplicate-value weights, missing/infinite data and short/empty inputs
- 122 independent checks: 59 broad invariants/edge cases plus 63 ordering/weight
  invariants, including explicit regressions for documented upstream defects

The one example exercises row means, column medians, weighted summaries,
quantiles, shape-preserving column ranks, set averages and a sum identity.

## Independent semantic checks

A separately compiled bounds-checked probe tested 1,500 deterministic generated
weighted samples through ten scalar summaries, for 15,000 comparisons. It found
130 upstream differences, all explained by the documented odd-count
infinite-weight median defect; the mathematically corrected oracle matched all
15,000 outputs. There were no other weighted numerical mismatches.

The documented omitted-weight subset defects were also independently reproduced
against the pinned R package. These deliberately repaired cases are not counted
as exact matches to upstream. See README.md, API_COVERAGE.md and per-function FPM
notes for the exact examples and limitations. Extreme finite overflow/cancellation
improvements are likewise explicit deviations, not claims of bitwise compatibility.

## Reproducible fixtures

Both R fixture generators were rerun separately and matched the committed
text files byte for byte:

| Fixture | SHA-256 |
|---|---|
| test/reference_values.txt | ea491a0353d87a918081d9f9c696148de2032e67373c3fb545de86a2a137bc3a |
| test/order_reference.txt | 411312e2aacf597525b7adaf33e6f7fb117538d8b13eb82d75ac5829fbe6b65b |

`tools/generate_references.py` records every R expression, checks the package
version and regenerates the embedded end-to-end Fortran oracle. Running
`Rscript test/generate_order_reference.R` regenerates the weighted/order fixture.
Ordinary builds and tests need neither R nor Python.

## Source and coverage audit

Static audit passed with zero errors and zero warnings:

- 99 function mapping entries, 99 computational exports, fraction 1.0 and empty
  untranslated list; README and manifest agree
- Three allocation-only R exports explicitly excluded from the denominator
- Every mapped public Fortran procedure and recorded source path exists
- All Fortran is free-form; one `dp` real kind is imported from r_kinds
- Every dummy has explicit INTENT or VALUE, its own declaration and meaningful
  trailing FORD documentation
- No duplicate maintained Fortran source, prohibited semicolon-separated
  statements, self-comparison NaN tests, vendored shared dependencies or system
  BLAS/LAPACK linkage
- NaN/Inf behavior is tested under ordinary debug and release optimization
- Complete Artistic-2.0 and GPL-2.0 texts and preserved provenance/attribution

Before packaging, run `fpm clean --all` and ensure the archive has no build
products, module/object files, executables, caches or nested archives. The
package directory must remain named exactly `matrixStats`.
