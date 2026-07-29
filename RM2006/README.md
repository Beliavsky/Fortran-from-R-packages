# RM2006 Fortran

A modern Fortran translation of the CRAN package **RM2006 0.1.1**, which
implements the RiskMetrics 2006 conditional covariance methodology described
by Zumbach (2007).

The library takes a `T x K` matrix of returns and produces a `K x K x (T+1)`
array. The first slice is the multiscale backcast. Slice `t+1` is the covariance
estimate after incorporating observation `t`, so the last slice is the
one-step-ahead covariance forecast after the complete sample.

## Features

- Full multiscale RiskMetrics 2006 recursion.
- Original defaults: `tau0=1560`, `tau1=4`, `kmax=14`, `rho=sqrt(2)`.
- Public `rm2006` compatibility alias and descriptive
  `rm2006_covariance` procedure.
- Public scale and weight calculation through `rm2006_scale_weights`.
- Input validation and status codes instead of abrupt termination.
- No external numerical-library dependencies.
- Deterministic numerical tests against an independent implementation.
- Safe behavior when the sample contains fewer observations than `kmax`.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example custom_parameters
```

The package uses standard Fortran 2018 and has no dependencies beyond an FPM
compatible Fortran compiler.

## Minimal example

```fortran
program example_rm2006
   use rm2006_kinds, only : dp
   use rm2006_module, only : rm2006
   implicit none

   real(dp) :: returns(100, 3)
   real(dp), allocatable :: covariance(:, :, :)

   ! Fill returns(t, asset) here.
   returns = 0.0_dp

   call rm2006(returns, covariance)
   print *, covariance(:, :, size(covariance, 3))
end program example_rm2006
```

Optional arguments may be supplied by keyword:

```fortran
call rm2006(returns, covariance, tau0=1560.0_dp, tau1=4.0_dp, &
            kmax=14, rho=sqrt(2.0_dp), status=status)
```

See `API.md`, `PORTING.md`, and `TESTING.md` for further details.

## License

The original R package declares `GPL (>= 2)`. This translation is therefore
licensed under **GPL-2.0-or-later**. Complete GPL-2.0 and GPL-3.0 texts are
included. Original source and metadata are retained under `original/` for
provenance.

## Reference

Gilles Zumbach (2007), *The RiskMetrics 2006 Methodology*, SSRN 1420185.
