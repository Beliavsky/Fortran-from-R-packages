# uncorbets-fortran

Modern Fortran implementation of the computational routines in the R package
`uncorbets`: minimum-torsion diversification, PCA torsion, effective-number-
of-bets analysis, and ENB-maximizing portfolio weights.

## Features

- Principal-components torsion
- Approximate and exact minimum-torsion matrices
- Symmetric positive-semidefinite matrix square roots
- Diversification probabilities and effective number of bets
- Long-only, fully invested ENB maximization
- Analytical ENB gradient and numerical Hessian diagnostics
- Typed status and result objects
- Self-contained linear algebra and optimization
- FPM and Makefile builds

## FPM

```text
fpm build
fpm test
fpm run
fpm run --example basic_uncorbets
```

## Make

```text
make check
make optimized
```

## Minimal example

```fortran
program example
  use uncorbets, only : dp, torsion_result, effective_bets_result, &
      torsion, effective_bets
  implicit none

  real(dp) :: sigma(3, 3), weights(3)
  type(torsion_result) :: tresult
  type(effective_bets_result) :: bets

  sigma = reshape([1.00_dp, 0.35_dp, 0.15_dp, &
                   0.35_dp, 1.50_dp, 0.25_dp, &
                   0.15_dp, 0.25_dp, 2.00_dp], [3, 3])
  weights = 1.0_dp / 3.0_dp

  tresult = torsion(sigma, model='minimum-torsion', method='exact')
  bets = effective_bets(weights, sigma, tresult%matrix)
  write(*, *) bets%probability
  write(*, *) bets%enb
end program example
```

## License

MIT, matching the upstream package. The upstream R package calls the
GPL-licensed `NlcOptim` package for one constrained optimization. This port
uses an independently implemented projected-gradient simplex optimizer, so no
GPL code is compiled or linked and the MIT license is retained.
