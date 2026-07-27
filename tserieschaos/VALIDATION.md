# Validation

## Environment

- Compiler: GNU Fortran 14.2.0
- Language mode: Fortran 2018
- External numerical libraries: none
- Operating environment: Debian Linux container

## Strict debug build

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

## Optimized build

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror
-fbacktrace
```

## Test coverage

### Core embedding and nonlinear statistics

- Regular delay embedding with `m=3`, `d=2`
- Explicit lag embedding
- Multivariate embedding and column ordering
- Hand-derived correlation integral equal to `2/3`
- Multi-dimension correlation curves with monotonicity and terminal-value checks
- Lag-zero AMI equal to `log(2)` for an alternating binary sequence
- Finite-sample lag-one AMI using the original package's marginal convention
- Exact normalized recurrence-distance matrix
- Space-time isoline dimensions, non-negativity, and quantile ordering

### Direct, box-index, and automatic neighbor search

- Search-method selection and reporting for `direct`, `box`, and `auto`
- Exact equality of direct and box k-nearest indices
- Direct and box distance equality to `1e-14`
- Deterministic ordering by distance and sample index
- Self and Theiler-window exclusion
- Original-scale returned distances
- Exact equality of direct and box false-neighbor counts and fractions
- Exact equality of direct and box Lorenz Lyapunov reference counts
- Direct and box Lyapunov stretching-path agreement to `1e-14`
- Automatic box selection for sufficiently large low-dimensional embeddings
- Invalid and unavailable neighbor slots remain represented by `-1`

A deterministic 1,200-point performance-regression case records exact distance calculations rather than unstable wall-clock timing:

```text
k-nearest direct evaluations:       474415
k-nearest box evaluations:           87901
false-neighbor direct evaluations:  1410156
false-neighbor box evaluations:      261242
```

The box method therefore performs fewer than one quarter of the direct exact-distance evaluations in both tested workflows while returning identical numerical results.

### Other neighbor and Lyapunov routines

- False-neighbor curves over multiple embedding dimensions
- k-nearest search shape and distance ordering
- Exact zero log-stretching for two trajectories separated by a constant distance
- Exact linear-fit intercept `2` and slope `0.5`
- End-to-end Lyapunov neighbor tracking on a Lorenz trajectory

### Continuous systems

- Hand-calculated Lorenz, Rossler, and Duffing derivatives
- Lorenz equilibrium preservation under RK4
- Duffing zero orbit with exact linear phase evolution
- Finite Rossler trajectory
- Custom observation callback applied to an integrated trajectory

### Executables

The following are compiled and run in both build modes:

- `demo_tserieschaos`
- `analyze_csv data/lorenz_sample.csv`
- `analyze_csv data/lorenz_sample.csv box`
- `analyze_csv data/lorenz_sample.csv direct`
- `lorenz_analysis`

### Licensing

`test/check_license.sh` verifies:

- the complete GPL version 2 license file,
- `GPL-2.0-only` in `fpm.toml`, and
- SPDX and GPL-2-only notices in every `.f90` file under `src`, `app`, `example`, and `test`.

## Results

```text
Core embedding and nonlinear-statistic tests passed.
Direct, box-index, auto-search, and Lyapunov tests passed.
Continuous-system and RK4 tests passed.
GPL-2.0-only source license checks passed.
debug build, tests, and applications passed.
release build, tests, and applications passed.
```

## Claims not made

- Exact equivalence to `deSolve::lsoda`
- Exact R random streams or R time-series metadata
- Identical memory layout or box collision behavior to the original fixed C arrays
- Identical floating-point tie handling at `d2` bin boundaries
- Plotting or S3 behavior
- Validation through `fpm`, because `fpm` was unavailable
