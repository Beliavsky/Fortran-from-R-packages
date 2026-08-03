# PINstimation modern Fortran

A self-contained modern Fortran/FPM translation of the computational core of
[`PINstimation`](https://github.com/monty-se/PINstimation) 0.2.0.

The library estimates the probability of informed trading under the original
PIN model, multilayer PIN (MPIN), adjusted PIN (AdjPIN), and volume-synchronized
PIN (VPIN/IVPIN). It also provides simulation, posterior probabilities, layer
detection, and high-frequency trade classification.

## Main capabilities

- Stable direct, E, LK, and EHO PIN log-likelihoods
- PIN MLE, adapted EA/GWJ/YZ initialization paths, posterior probabilities,
  good/bad decomposition, simulation, and random-walk Bayesian sampling
- MPIN likelihood, MLE, ECM estimation, state posteriors, simulation, and layer
  detection
- AdjPIN six-state likelihood, restricted/unrestricted MLE and ECM, simulation,
  AdjPIN and PSOS measures
- Tick, quote, Lee-Ready, and Ellis-Michaely-O'Hara classification
- Volume-bucket construction, VPIN rolling measures, and rolling IVPIN MLE
- Self-contained Poisson generation, numerical optimization, and utilities

All public numerical inputs use `real(dp)`, where `dp = kind(1.0d0)`. Daily
buy/sell counts use 64-bit integers.

## Build with FPM

```text
fpm build
fpm test
fpm run --example pin_estimation
fpm run demo_pinstimation
```

The manifest version is numeric and FPM-compatible: `0.2.0`.

## Build directly with GNU Fortran

On Unix-like systems:

```text
./scripts/build_all.sh check
./scripts/build_all.sh release
```

On Windows with `gfortran` available in `PATH`:

```text
scripts\build_all.bat check
scripts\build_all.bat release
```

## Minimal example

```fortran
use pinstimation

type(pin_parameters) :: truth
type(pin_result) :: estimate
type(trade_counts) :: data
integer, allocatable :: state(:)

truth = pin_parameters(0.35_dp, 0.45_dp, 15.0_dp, 20.0_dp, 18.0_dp)
call simulate_pin(300, truth, data, state, seed=2026)
call fit_pin(data, estimate)
print *, estimate%pin
```

## Translation boundaries

The native API uses arrays and derived types rather than R formulas, data
frames, S4 classes, names, progress bars, futures, or tidyverse pipelines. The
VPIN path accepts already aggregated time-bar price changes, volumes, and
durations; calendar parsing and insertion of empty bars are intentionally left
to the application. See `docs/API_MAP.md` and `docs/PORTING_NOTES.md` for exact
and adapted coverage.

## License and provenance

The upstream package is licensed under GPL version 3 or later. This translation
uses the same license. The full upstream source snapshot is retained under
`upstream/PINstimation-master`.
