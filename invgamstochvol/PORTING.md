# Porting notes

## Directly translated numerical methods

The following computations follow `src/invgamstochvol.cpp`:

- logarithms of rising factorials;
- truncated `2F1` hypergeometric sums;
- the observation-by-observation exact likelihood recursion;
- the coefficient and scale tables returned for smoothing;
- backward mixture-index sampling and conditional gamma draws.

A fixed eight-observation regression case agrees with the original C++ source
to approximately `2e-15` in the total likelihood and each likelihood
contribution.

## Numerical changes

The Fortran implementation uses log-sum-exp reductions when combining positive
series terms. This is algebraically equivalent to the upstream maximum-plus-7
scaling but is less likely to overflow or underflow.

The direct hypergeometric recurrence handles `zstar = 0` explicitly. The C++
implementation reaches the same mathematical value through logarithms of zero.

## Random numbers

The original package seeds a C++ Mersenne Twister from the system clock for
each call. The Fortran API accepts an optional 64-bit seed and uses a portable
Park-Miller uniform stream, Box-Muller normals, and Marsaglia-Tsang gamma
sampling. Repeated calls with the same seed are reproducible across supported
compilers. Omitting the seed uses a documented fixed default, which favors
reproducible scientific tests over clock-dependent behavior.

## Parallel arguments

`nproc` and `nproc2` remain optional arguments for API familiarity. The current
FPM library is serial and dependency-free; positive values are accepted and do
not alter results.

## Data and plotting

The bundled R data object is retained under `original/`. The package contains
no computational plotting API. Plotting shown in the R examples is intentionally
left to the calling application.
