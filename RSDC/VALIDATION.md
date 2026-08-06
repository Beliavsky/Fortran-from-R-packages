# Validation

The checked suite is compiled with:

```text
-std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace -O0
```

The optimized suite is compiled with:

```text
-std=f2018 -Wall -Wextra -Werror -Wno-maybe-uninitialized -O3
```

Tests cover:

1. canonical partial-correlation round trips and transition maps;
2. fixed-P and TVTP Hamilton filtering, smoothing, normalization, and likelihood parity;
3. simulation, covariance/correlation forecasts, multi-step forecasts, and Viterbi decoding;
4. constant, fixed-transition, and TVTP estimation with regime ordering;
5. minimum-variance and maximum-diversification portfolios, warm starts, and ergodic diagnostics;
6. numerical scores, OPG covariance, correlation bands, and parametric bootstrap.

`TEST_RESULTS.txt` records the final checked and optimized runs.

The optimized build disables only `-Wmaybe-uninitialized`, because GCC 14 emits
known false positives for unallocated derived-type descriptors in the vendored
DEoptimR optional-constraint branches. Runtime-checked builds retain all warnings.
