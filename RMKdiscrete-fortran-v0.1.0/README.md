# RMKdiscrete-fortran

Modern Fortran/FPM translation of the computational code in the R package
`RMKdiscrete` 0.1 by Robert M. Kirkpatrick.

The upstream package describes itself as providing "Sundry discrete probability
distributions and helper functions" and is licensed GPL (>= 2). This translation
is therefore distributed as GPL-2.0-or-later. See `LICENSE`, `COPYING`, and
`doc/UPSTREAM_DESCRIPTION`.

## Implemented functionality

The public `rmkdiscrete` module provides:

- Lagrangian Poisson/generalized-Poisson distribution:
  - `dlgp`, `plgp`, `qlgp`, `rlgp`, `rlgp_sample`
  - `lgp_findmax`, `lgp_get_nc`, `slgp`
  - explicit equivalents of all valid `LGPMVP` parameter conversions
- Negative-binomial helpers:
  - `dnegbin`, `rnegbin`, `rnegbin_sample`
  - explicit equivalents of all valid `negbinMVP` conversions
- Bivariate Lagrangian Poisson:
  - `dbilgp`, `rbilgp`, `rbilgp_sample`, `bilgp_logmv`
- Bivariate negative binomial:
  - `dbinegbin`, `rbinegbin`, `rbinegbin_sample`, `binegbin_logmv`
- Mana Clash distributions:
  - `dmanaclash_xyn`, `dmanaclash_dmg`, `dmanaclash_net`
  - `rmanaclash`, `rmanaclash_sample`

The bivariate distributions retain the upstream shared-component construction:
`Y1 = U + X` and `Y2 = U + Y` with independent latent counts.

## Build with FPM

```text
fpm build
fpm test
fpm run --example demo_rmkdiscrete
```

No external Fortran dependencies are required.

## Small example

```fortran
program example
  use rmkdiscrete, only : dp, dlgp, plgp, qlgp
  implicit none

  print *, dlgp(3, 2.0_dp, 0.3_dp)
  print *, plgp(3.0_dp, 2.0_dp, 0.3_dp)
  print *, qlgp(0.5_dp, 2.0_dp, 0.3_dp)
end program example
```

## Translation notes

See `TRANSLATION_NOTES.md` for the R-to-Fortran API map, numerical details, and
intentional implementation differences.
