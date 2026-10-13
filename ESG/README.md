# ESG: economic scenarios and asset distributions

Modern free-form Fortran/FPM translation of CRAN ESG 1.3. It implements
Hull–White short-rate paths, correlated lognormal stocks, real estate, liquidity
spreads, the upstream default-spread recurrence, martingale diagnostics, and nine
asset pricing distributions. All computational procedures are pure (including
caller-owned RNG state advancement); scalar model kernels are pure elemental.

**This is a source-compatible research translation, not a corrected or validated
financial model.** Several upstream formulas and time-index conventions are
unusual or mathematically questionable. They are retained and tested explicitly.
Do not infer production pricing suitability from passing regression tests.

## Build and use

Place this directory at repository root beside the existing `rfortran-core` and
`randompack` packages, with randompack's `rfortran-linalg` sibling present:

```
fpm build
fpm test
fpm run --example scenarios
```

A recent gfortran and FPM are required. The package uses portable free-form
Fortran, with no R, C, Python, system BLAS or system LAPACK runtime requirement.
The sibling linear-algebra package fetches its pinned Fortran LAPACK dependency
through FPM. First-time resolution needs access to its configured Git repository
(or a verified local FPM cache). No platform-specific link flags are introduced;
Linux gfortran execution has been tested, but this is not a claim of a Windows
runtime test. Do not compile with `-ffast-math`, `-Ofast`, `-ffinite-math-only` or
similar flags: documented quiet-NaN and nonfinite validation require IEEE
semantics. Dependency packages are not bundled in the ESG archive.

```fortran
use esg
implicit none
type(esg_parameters) :: p
type(esg_scenario) :: s
type(randompack_rng_type) :: rng
real(dp) :: zero(73), z(1000,4,5)
real(dp), allocatable :: prices(:)
integer :: status
zero = 0.03_dp
rng = randompack_rng('pcg64', seed=12345)
call esg_innovations(rng, z, status)
call esg_generate(p, zero, z, s, status)
if (status /= 0) error stop 'scenario generation failed'
call esg_zc_bond(s, 1, 3.5_dp, prices, status)
if (status /= 0) error stop 'pricing failed'
```

For real work, check status immediately after every call. `example/scenarios.f90`
is a complete program, including the innovation-status check.

## Model and array conventions

- `esg_parameters` contains named risk parameters; inspect the type declarations
  for defaults. They are example values, not calibrated market estimates.
- Innovations have shape `(number_of_scenarios, horizon, 5)`. The last dimension
  is short-rate, stock, estate, liquidity, default, in that order. Inputs may be
  supplied explicitly for controlled experiments and comparisons with R.
- Every path has shape `(number_of_scenarios, horizon+1)`. Column 1 is time zero;
  transitions are annual. The source comments sometimes call `t` a month, but
  the computations mix annual paths and monthly curve indices. This package
  does not reinterpret or harmonize those units.
- `esg_forward_rates` computes annual samples of one-month forward differences:
  sample i is `(12*i+1)*zero(12*i+1) - 12*i*zero(12*i)` for i=1..horizon+1.
  The monthly input therefore needs at least `12*(horizon+1)+1` rates.
  `esg_zero_rates` retains the first `12*horizon+1` monthly inputs.
- `esg_generate` defaults to all five risks; optional kind is `shortRate`,
  `stock`, `realEstate`, `liquiditySpread` or `defaultSpread`. Curves are required
  by this typed scenario constructor even for spread-only kinds; the standalone
  spread procedures need no curves. Unrequested paths remain unallocated.
- Standalone `esg_short_rate` accepts extracted annual forwards. `esg_stock`
  and `esg_real_estate` accept rate paths. No R S4 objects, recycled arguments,
  implicit global state, R warnings, or plotting are implemented.
- Valuation results have one element at t=0, otherwise one per scenario. Unlike
  R's scalar early returns at maturity, Fortran retains that result shape.
  `esg_asset_price` dispatches on the nine original asset type strings;
  `esg_simulate_asset` constructs the relevant risks and then values them.
- Volatilities are nonnegative, `k>0`, `|rho|<=1`, initial stock/estate positive,
  initial default spread/alpha/beta nonnegative. LMN corporate/CDS pricing
  requires strictly positive default volatility. Coupons must be positive
  integers, at most 10000, with at most one million payments; coupon payments
  are nonnegative, strikes positive, and omega in [0,1]. These explicit guards
  replace R's often incomplete or permissive validation.

Status values: 0 success; 1 invalid domain, nonfinite input or missing required
argument; 2 inconsistent/missing arrays or curve/horizon range; 3 a negative
default-spread path; 4 nonfinite/overflowed numerical output. Do not use outputs
when status is nonzero. Default generation retains the first negative column
and leaves following columns as quiet NaN, instead of evaluating sqrt(negative).
No positivity clipping is applied. Stock overflow gives quiet NaN/status 4;
exponentials below the normal floating-point range return zero. Memory allocation
failure is handled by the Fortran runtime rather than a package status code.

## Deliberately retained upstream quirks

1. Short-rate drift includes `vol^2/2 * K(i-1)^2 * (1-exp(-k))`. Stock and estate
   growth use the newly generated rate in column i+1, not the preceding rate.
2. Default drift subtracts `volDefault/4`, **not** its square divided by four.
   This can create negative intensities even with nonnegative inputs. Its drift
   is alpha*(beta-lambda), whereas the separate LMN pricing coefficients use
   alpha and beta in the original A/B/G/H formulas; the translation does not
   recalibrate or reconcile those different parameter conventions.
3. Stock option valuation at positive t uses stock column t (one period lagged),
   and Tau uses absolute maturity T rather than T-t. The unusual positive second
   term in upstream Tau is preserved. Stock options return one at t=T instead
   of intrinsic payoff. Convertible bonds add a stock call with strike one.
4. Zero-coupon prices use monthly rate index `12*T` at grid maturities, but
   discount the valuation date using rate index `t+1`. Off-grid prices interpolate
   **discount factors**, not rates. R's exact monthly-grid comparisons and
   truncating numeric subscripts are retained, including boundary-roundoff quirks.
   The first interpolated endpoint at zero is one.
5. Bond options use `log(Pts/(PtT*Strike)+H/2)/H`, with H/2 **inside** the logarithm,
   and H depends on absolute T and s. This is not the textbook Black formula.
6. Coupon frequency is not used to divide couponRate. The schedule includes a
   coupon on the valuation date when it lies on the coupon grid. Accrued interest
   contains the source's additional 3/360 term.
7. Corporate bonds retain bare `A(T)` for principal, rather than inserting the
   apparently omitted discount/survival/liquidity factors. Coupon LMN arguments
   use absolute payment dates. Corporate recovery is multiplied by `1-omega`;
   CDS uses `omega`. CDS does not use liquidity spreads in its numerical formula.
   Both retain 100 right-endpoint rectangular integration steps.

Numerical extensions: algebraically equivalent LMN A/B/G/H expressions avoid
upstream `exp(large)*tiny` overflow. Finite results beyond R's finite range are
an extension, not an R-verified financial correction. Small-argument K uses a
Taylor expansion. LMN coefficients use rationalized differences and core
log1p/expm1 helpers for tiny volatilities; a continuous zero-diffusion limit is
used only if a positive volatility squared underflows. A local Riccati-series
integral avoids cancellation when both beta and volatility are tiny. Zero option variance/H returns the continuous intrinsic-value
limit, where R may return NaN. The broken `rAssetDistribution` call spelling is
made usable: both `EuroCall_S` and `EuroCall_Stock` select the stock call. Source
models themselves are otherwise not repaired.

## Randomness and validation

`randompack_rng_type` provides independent caller-owned streams. `esg_innovations`
fills each risk block using that sibling's standard normals; fixed engine, seed,
array shape and call sequence are reproducible without changing Fortran's global
RNG. This is **algorithm parity with supplied shocks, not R RNG sequence parity**.
R uses a different generator. Its `rStock` and `rRealEstate` wrappers also generate
and discard an initial block of rate normals before regenerating rates. The R
reference driver explicitly accounts for those discarded draws. This Fortran
API consumes only the supplied arrays and does not discard hidden draws.

Tests contain 469 genuine ESG 1.3 reference values from 71 cases, covering the
20 exported computational operations, all five risks, every pricing method and
the dispatcher at t=0,1,2, off-grid maturities, and terminal conventions. The
original package is loaded, its imported rnorm is temporarily replaced with a
controlled finite innovation queue, and the actual exported functions/methods
are invoked. Nonfinite R references are refused. Tests do not need R installed.
Reproduce fixtures (with ESG 1.3 installed) using:

```
Rscript reference/generate_reference.R reference/expected.tsv
python reference/embed_reference.py
fpm test
```

Independent tests also check normal moments and stream isolation, analytic
short-rate, log-stock and default-spread moments, requested innovation covariance, liquidity
variance, discounted-stock means, coupon cash flows, option parity, all nine
end-to-end asset branches, degenerate/invalid shapes and domains, tiny mean
reversion, deep tails, default negativity and overflow. Fixed local streams make
statistical tests deterministic; tolerances allow expected finite sample error.

## Translation coverage

Status: **substantial**. **20 of 20 (100%)** exported computational R functions
are mapped. The fraction measures mapped computational functions, not complete
R compatibility, corrected finance theory, or identical stochastic streams.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| setForwardRates | esg_models | esg_forward_rates | substantial |
| rShortRate | esg_models | esg_short_rate | substantial |
| rStock | esg_models | esg_stock | substantial |
| rRealEstate | esg_models | esg_real_estate | substantial |
| rLiquiditySpread | esg_models | esg_liquidity | substantial |
| rDefaultSpread | esg_models | esg_default | substantial |
| rAllRisksFactors | esg_models | esg_generate | substantial |
| rAssetDistribution | esg_pricing | esg_simulate_asset | substantial |
| MartingaleTest | esg_models | esg_martingale | substantial |
| customPathsGeneration | esg_models | esg_generate | substantial |
| Asset_PriceDistribution | esg_pricing | esg_asset_price | substantial |
| Bond_PriceDistribution | esg_pricing | esg_bond | substantial |
| CDSPremium_PriceDistribution | esg_pricing | esg_cds | substantial |
| ConvBond_PriceDistribution | esg_pricing | esg_convertible_bond | substantial |
| CBond_PriceDistribution | esg_pricing | esg_corporate_bond | substantial |
| EuroCall_Stock_PriceDistribution | esg_pricing | esg_stock_call | substantial |
| EuroPut_Stock_PriceDistribution | esg_pricing | esg_stock_put | substantial |
| EuroCall_ZC_PriceDistribution | esg_pricing | esg_zc_call | substantial |
| EuroPut_ZC_PriceDistribution | esg_pricing | esg_zc_put | substantial |
| ZCBond_PriceDistribution | esg_pricing | esg_zc_bond | substantial |

See [API-COVERAGE.md](API-COVERAGE.md) for the exact denominator and exclusions,
and [NOTICE.md](NOTICE.md) for provenance and license attribution.
