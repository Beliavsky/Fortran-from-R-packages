# intradayModel-fortran

A self-contained modern Fortran translation of the computational core of R package
`intradayModel` 0.0.1.

The library fits and applies the univariate intraday-volume state-space model of
Chen, Feng, and Palomar. On the log scale, the observed signal is decomposed into:

- a daily latent component that is constant within each trading day;
- a dynamic intraday AR(1) component;
- a deterministic periodic seasonal profile; and
- measurement noise.

## Implemented API

- `fit_volume`: ordinary EM or accelerated EM fitting, with any parameter fixed or initialized by the caller.
- `decompose_volume`: Kalman-smoothed analysis or one-bin-ahead forecasting.
- `forecast_volume`: convenience wrapper for forecast decomposition.
- `uniss_kalman` and `uniss_em_update`: lower-level filter, smoother, and EM routines.
- `simulate_intraday_volume`: deterministic portable simulation helper for tests and examples.

The main model and output structures are `volume_parameters`, `volume_model_spec`,
`volume_fit_control`, `volume_model`, and `volume_decomposition`.

## Build with FPM

```text
fpm build
fpm test
fpm run --example basic_fit
fpm run --example forecast_example
fpm run demo_intraday_model
```

The project has no external numerical dependencies.

## Direct GNU Fortran validation

On Unix-like systems:

```text
./scripts/test_gfortran.sh
./scripts/test_gfortran_optimized.sh
```

On Windows with GNU Fortran available in `PATH`:

```text
scripts\test_gfortran.bat
```

## Data layout

Input volume is a positive real matrix with shape `n_bin x n_day`. Fortran column-major
flattening therefore follows chronological order: all bins of day 1, then all bins of
day 2, and so on. Columns containing a nonfinite or nonpositive value are removed by
`clean_volume_data`.

## Scope

The compiled library omits `generate_plots`, `ggplot2`, `patchwork`, `xts`, and `zoo`
infrastructure. Dates and time-zone-aware indexing remain the caller's responsibility.
The original package sources are retained under `original/`.

## License

Apache License 2.0, matching the original package metadata. See `LICENSE` and
`NOTICE.md`.
