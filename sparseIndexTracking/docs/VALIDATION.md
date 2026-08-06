# Validation

The package is compiled with GNU Fortran 14.2 using two configurations.

Checked configuration:

```text
-std=f2018 -Wall -Wextra -Werror -O0 -g -fcheck=all -fbacktrace
```

Optimized configuration:

```text
-std=f2018 -Wall -Wextra -Werror -O3
```

## Tests

### `test_projection`

Checks a capped-simplex projection against an analytic four-asset solution,
verifies the direct `bisection` interface, constraints, and infeasible caps.

### `test_ete_recovery`

Constructs a benchmark exactly from two of six deterministic assets. The ETE
solver recovers weights close to 0.7 and 0.3, removes the other assets, and
matches the `spIndexTrack` compatibility entry point. It also checks the
upstream threshold-and-renormalize behavior.

### `test_tracking_measures`

Fits ETE, DR, HETE, and HDR to contaminated deterministic returns. Every result
is feasible, and HETE is closer than ETE to the uncontaminated generating
portfolio after a large benchmark outlier.

### `test_objectives`

Verifies exact hand-computed ETE, DR, HETE, corrected HDR, and
source-compatible HDR objective values.

### `test_errors`

Checks univariate input, zero lambda, infeasible caps, missing Huber parameters,
unknown measures, invalid initial portfolios, and safe iteration-limit return.

## Example

The example benchmark is generated from three of eight assets and contaminated
with one large outlier. The HETE fit selects the three generating assets and
reports a finite tracking RMSE.

An R executable was not available in the build environment, so validation is
based on formula-level fixtures, analytic projection results, deterministic
recovery cases, and strict checked/optimized Fortran builds rather than a live
R-to-Fortran numerical comparison.
