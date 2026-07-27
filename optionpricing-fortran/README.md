# optionpricing-fortran

Modern Fortran translation of the numerical routines in the R package
`OptionPricing` 0.1.2 by Wolfgang Hormann and Kemal Dingec.

The library prices European and discretely monitored arithmetic-average Asian
options under geometric Brownian motion. It includes analytic control-variate
formulas, Lord's approximation, ordinary Monte Carlo, conditional Monte Carlo,
multiple control variates, randomized Korobov lattice rules, and estimators of
price, delta, and gamma.

## Build

```text
fpm build
fpm test
fpm run optionpricing_demo
fpm run --example european_and_analytic
fpm run --example monte_carlo_methods
```

No external numerical or statistics library is required.

## Main API

```fortran
use optionpricing

type(european_result) :: euro
type(greeks_result)   :: asian

euro = bs_ec(0.25_dp, 100.0_dp, 0.05_dp, 0.2_dp, 100.0_dp)

asian = asian_call(1.0_dp, 12, 100.0_dp, 0.05_dp, 0.2_dp, 100.0_dp, &
   "best", "QMC", n=2039, nout=50, a_gen=1487, seed=4711)
```

`greeks_result%estimate` and `greeks_result%error95` are ordered as price,
delta, and gamma.

Important public procedures include:

- `bs_ec`, `bs_ep`, `bs_european_call`, `bs_european_put`
- `asian_call_app_lord`
- `eval_ecv`, `eval_lb`, and `eval_eqcv`
- `asian_call_naive_mc`, `asian_call_ncv_lr_mc`, `asian_call_cmc_cv`
- `asian_call_best_mc`, `asian_call_naive_qmc`, `asian_call_best_qmc`
- `korobov_lattice`, `naive_pca_matrix`, and
  `conditional_generation_matrix`
- `conditional_estimates_z` for caller-supplied normal points

## Numerical design

- Double precision is defined as `dp = kind(1.0d0)`.
- Linear regressions use rank-revealing, column-pivoted QR rather than normal
  equations.
- Symmetric eigendecomposition is self-contained.
- Normal CDF, density, and inverse CDF routines are self-contained.
- Randomized QMC uses the Korobov lattice supported by the original package.
- Optional integer seeds initialize Fortran's intrinsic random generator.

The upstream Sobol branch only printed that Sobol was not activated and its
actual call was commented out. It is therefore not represented as a working
method in this port.

## Compatibility correction

The original `BS_EC` and `BS_EP` functions calculate gamma with `pnorm(d1)`.
Black-Scholes gamma requires `dnorm(d1)`. The Fortran result field `gamma`
contains the correct value, while `upstream_gamma` preserves the original R
calculation for comparison.

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for details.

## License

The supplied package declares `GPL-2 | GPL-3`. This project preserves that
choice as `GPL-2.0-only OR GPL-3.0-only`. Complete copies of GPL version 2 and
GPL version 3 are included. The original source is retained unmodified under
`original/OptionPricing-0.1.2`.
