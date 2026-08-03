# JumpTest modern Fortran translation

This project translates the computational code of the R package **JumpTest 1.1**
to standard-conforming modern Fortran with an FPM project layout.

## Implemented routines

### Stochastic-volatility simulation

- `sv`: one-factor square-root stochastic volatility without jumps.
- `svj`: square-root stochastic volatility with compound Poisson jumps.
- `sv1f`: Chernov-style one-factor log-volatility model.
- `sv1fj`: one-factor log-volatility model with jumps.
- `sv2f`: two-factor log-volatility model.

Every simulator accepts an optional 64-bit seed and returns a
`simulation_result` containing the price path, latent factor path, jump counts,
jump sizes, and status information where applicable.

### Jump tests and p-value pooling

- `bns_statistic`, `amin_statistic`, and `amed_statistic`.
- `jumptestday` and `jumptestperiod`.
- `pcombine` for assembling method-specific p-value matrices.
- `ppool` with `SD`, `FD`, `SI`, `FI`, `MI`, and `MA` pooling.
- Benjamini-Hochberg adjustment through `bh_adjust`.

## Build with FPM

```text
fpm build
fpm test
fpm run demo_jumptest
```

Direct GNU Fortran validation scripts are provided in `scripts/`.

## Minimal example

```fortran
use jumptest, only : dp, statp_result, jumptestday

real(dp) :: returns(5)
type(statp_result) :: result

returns = [0.001_dp, -0.002_dp, 0.012_dp, -0.001_dp, 0.0005_dp]
call jumptestday(returns, result, 'BNS')
print *, result%stat, result%pvalue
```

## Interface notes

Fortran is case-insensitive, so the R arguments `M` and `m` cannot coexist as
separate names. Simulator interfaces use `intervals_per_period` and `periods`.
The total number of simulated observations is their product, and the time step
is `1 / intervals_per_period`.

The upstream `SV` function returns the initial price, whereas its other four
simulators remove it. This inconsistency is preserved for compatibility.

The upstream `SVJ` and `SV1FJ` wrappers use malformed `sapply(..., rnorm, n=1,
mean=0)` calls. The Fortran port implements the apparent intended models:
`SVJ` uses a jump standard deviation of `sqrt(count*sigma1)`, and `SV1FJ`
uses the Poisson count as the jump standard deviation. See `PORTING.md`.

## License

The upstream package uses the MIT license. This translation is also distributed
under the **MIT license**. Original package sources and metadata are retained
under `original/JumpTest/`.
