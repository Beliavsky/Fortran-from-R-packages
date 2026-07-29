# Independent Reference Generation

Fixed distribution references were generated independently in Python using
SciPy special functions and root solving. The reference calculations use the
published Burr, generalized-gamma, generalized-F, q-Weibull, mixture, and
Birnbaum-Saunders formulas rather than calling the Fortran implementation.

Model recursion tests use hand-computed conditional means. Transaction-duration
references use explicitly enumerated timestamps and events. Diurnal tests include
synthetic Fourier curves for which ordinary least-squares recovery is exact up to
floating-point rounding.

The ACD estimation test simulates from known ACD(1,1) parameters with the
project's explicit RNG and verifies parameter recovery, likelihood outputs,
inference matrices, and forecasts. It is a recovery test rather than a claim of
bit-for-bit agreement with R's optimizer.
