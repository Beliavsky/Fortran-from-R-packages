# Jdmbs modern Fortran translation

This project translates the computational code of the R package **Jdmbs 1.4**
to standard-conforming modern Fortran with an FPM project layout.

## Implemented routines

- `normal_bs`: Monte Carlo geometric-Brownian option pricing.
- `jdm_bs`: independent jump-diffusion simulation for each asset.
- `jdm_new_bs`: shared jump events transmitted across companies through a matrix.

Every routine returns call and put estimates, Monte Carlo standard errors,
terminal prices, status information, the event count, and optionally complete
paths.

## Upstream compatibility and corrected scaling

The original R implementation draws one standard normal per day but does not
multiply it by `sqrt(dt)`. It also uses an off-by-one time index in
`normal_bs` and `jdm_bs`. The Fortran control type therefore provides:

```fortran
control%legacy_mode = .true.   ! default: reproduce upstream formulas
control%legacy_mode = .false.  ! standard daily GBM scaling
```

Corrected mode evolves prices recursively using

```text
S(t+dt) = S(t) exp((mu - sigma^2/2) dt + sigma sqrt(dt) Z)
```

and applies every jump occurring in an interval. Legacy mode preserves the
upstream absolute-price formulas and applies at most one jump per daily step.

The upstream package returns undiscounted expected payoffs. This remains the
default (`discount_rate = 0`), but a continuously compounded discount rate can
be supplied in `jdmbs_control`.

## Build with FPM

```text
fpm build
fpm test
fpm run demo_jdmbs
```

Direct GNU Fortran validation scripts are in `scripts/`.

## License

The upstream package declares `GPL (>= 2)`. This translation is distributed as
**GPL-2.0-or-later**. Original package sources and metadata are retained under
`original/Jdmbs/`.
