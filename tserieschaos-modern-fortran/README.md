# tseriesChaos Modern Fortran

A modern Fortran 2018 translation of the computational routines in the R package `tseriesChaos` 0.1-13.1.

The original package was written by Antonio Fabio Di Narzo and is licensed under GPL-2. This translation preserves that licensing as `GPL-2.0-only`.

## Implemented numerical features

- Univariate and multivariate delay embedding with regular or explicit lag sets
- Sample correlation integral at one scale
- Correlation-integral curves for multiple scales and embedding dimensions
- Histogram-based average mutual information
- False-nearest-neighbor fractions and neighbor counts
- Radius-limited k-nearest-neighbor searches with a Theiler window
- Direct exhaustive, box-index, and automatic neighbor-search modes
- Deterministic neighbor ranking by distance and then sample index
- Kantz-style neighbor following and logarithmic stretching paths
- Linear regression estimates from selected Lyapunov stretching intervals
- Normalized recurrence-distance matrices
- Space-time separation probability isolines from 10 percent through 100 percent
- Lorenz, Rossler, and Duffing differential systems
- Generic fourth-order Runge-Kutta integration
- Optional observation callbacks for simulated state trajectories

All routines return ordinary Fortran arrays, scalars, or status codes. No plotting or R object system is required.

## Neighbor-search modes

`false_nearest_fraction`, `false_nearest_curve`, `find_k_nearests`, and `lyapunov_stretching` accept the optional keyword:

```fortran
search_method="direct"  ! exhaustive reference implementation
search_method="box"     ! accelerated box-index implementation
search_method="auto"    ! select from data size, dimension, and radius
```

The default is `auto`. Automatic selection uses box search when there are at least 256 searchable points, the embedding dimension is at most 8, and the scaled radius is below 0.5. Otherwise it uses direct search.

The box index follows the strategy of the original C package: it bins points by the first and last coordinates of the embedding, scans neighboring boxes, and then applies the exact full-dimensional Euclidean-radius test. The Fortran version uses a dynamically populated grid of up to 100 by 100 cells instead of the original fixed modulo-indexed arrays. Therefore it avoids box collisions while preserving the neighbor set.

`find_k_nearests`, `false_nearest_fraction`, and `lyapunov_stretching` can optionally return `distance_evaluations` and `method_used`. The count is useful for deterministic performance regression tests and algorithm comparisons.

## Deliberate differences from R

- `sim.cont` in R delegates to `deSolve::lsoda`. The Fortran project provides a tested fixed-step RK4 integrator. It does not claim LSODA-equivalent trajectories or adaptive error control.
- The original C box search uses fixed 100 by 100 arrays and modulo binning. The Fortran box option uses an exact dynamic grid with at most 100 cells per coordinate. Direct and box paths are tested to return identical sorted neighbors and downstream statistics.
- The original `d2` code bins squared distances and then cumulatively sums the bins. The Fortran implementation directly counts each requested threshold. Results are mathematically equivalent away from exact bin boundaries; floating-point tie behavior can differ.
- Time-series dates, frequencies, windows, `ts` classes, S3 methods, and R's missing-value representation are not reproduced.
- Plotting and printing methods are excluded. The numerical matrices that those methods consume are included.

## Build and test

GNU Fortran:

```sh
make check
make check-opt
```

`make check` uses:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

`make check-opt` uses:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror
-fbacktrace
```

An `fpm.toml` manifest is included:

```sh
fpm test
fpm run demo_tserieschaos
fpm run analyze_csv -- data/lorenz_sample.csv auto
```

The manifest was not executed during release validation because `fpm` was unavailable in the validation environment.

## Applications

### Demonstration

```sh
./build/debug/demo_tserieschaos
```

Simulates the Lorenz system and computes AMI, a correlation integral, false-neighbor fractions, and a Lyapunov stretching slope. It reports the neighbor-search method selected by `auto`.

### CSV analysis

```sh
./build/debug/analyze_csv data/lorenz_sample.csv auto
./build/debug/analyze_csv data/lorenz_sample.csv box
./build/debug/analyze_csv data/lorenz_sample.csv direct
```

The application accepts either:

- one numeric value per row, or
- `Date,Value` rows with an optional header.

It reports sample scale, AMI, a correlation integral, false-neighbor fractions, the selected search method, and a Lyapunov estimate when the neighborhood contains enough points.

## Main modules

- `chaos_embedding`: delay-coordinate matrices
- `chaos_systems`: ODE systems and integration
- `chaos_metrics`: correlation sums, AMI, recurrence, space-time separation
- `chaos_neighbors`: direct/box neighbor search, false neighbors, and Lyapunov calculations
- `tserieschaos`: convenience umbrella module

See `API_MAP.md` for the mapping from each R/C routine to its Fortran counterpart and `VALIDATION.md` for the tested cases.
