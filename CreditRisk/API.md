# API and coverage

All 13 exported computational R functions are mapped in the manifest and
README. `cdsdata` is a dataset, not a computational export. The Fortran API
omits R object classes, recycling, interactive/plotting interfaces and the
dataset. Shapes and validity are explicit.

## Results and numerical domains

`use creditrisk` exports `dp`, `credit_curve_t`, `credit_ok=0`,
`credit_invalid=1`, and the pricing/simulation routines below.
`credit_curve_t%values` is an allocated matrix with dates down rows and the
following columns. `status` must equal `credit_ok` before the table is used.
For invalid inputs the result is invalid and is usually NaN-filled except
for dates; do not rely on other columns. Nonfinite output makes status invalid.
Structural curves require initial firm value above a strictly positive barrier.
Survival underflow to zero is treated as an invalid curve because log hazards
are no longer finite. No high-precision or logarithmic-output API is supplied.

| Procedure | Columns of `%values` |
|---|---|
| `merton` | maturity, firm value, equity value, debt value, survival |
| `black_cox`, `at1p` | maturity, firm value, barrier, survival, interval hazard |
| `sbtv` | maturity, survival, interval hazard |
| `cds`, `cds2` | maturity, survival, premium leg, protection leg, rate, price |
| simulation routines | One column per simulated path; no date column |

Every maturity or interval length is in years; rates, spreads and volatilities
are annual decimal values. Monetary inputs use consistent arbitrary currency
units. Default CDS premium is 0.005 and recovery is 0.4. Simulation allows a
first date of zero; other maturity vectors must be nonempty, strictly
increasing and positive. Negative finite interest rates are supported.
All array arguments must have their documented length. No R recycling occurs.

## Pricing and structural routines

- `cum_normal_density(x)`: pure elemental Hastings normal CDF approximation
  used by upstream option pricing. Infinite tails map to zero/one; NaN
  propagates. Tails beyond ±40 saturate before evaluating the polynomial.
  This is the source approximation, not exact normal probability; at zero
  it is about 0.500000000525, so exact identities involving that boundary
  can have the source's small approximation error.
- `generalized_black_scholes(option,spot,strike,time,rate,carry,sigma)`: pure
  elemental European call (`'c'`) or put (`'p'`) value. Spot/strike must be
  positive, time/volatility nonnegative, and all numeric inputs finite.
  Exact deterministic/expiry limits extend the R implementation. An invalid
  flag/domain or unrepresentable discounted asset/strike leg gives quiet NaN.
  Logarithm differences avoid unnecessary overflow in spot/strike ratios.
- `merton(debt,v0,sigma,rate,time)`: positive debt, firm value and volatility;
  scalar rate. This preserves the upstream survival numerator
  `log(v0/debt)+rate-sigma^2*time/2`, **not** the textbook numerator with
  `rate*time`. Equity and debt retain the upstream option calculation with
  expiry at the last time. Terminal at-the-money equity is explicitly zero.
- `black_cox(debt,v0,sigma,rate,gamma,time,barrier)`: optional `barrier` is K,
  default debt, with `0<K<=debt`. Initial barrier is `K*exp(-gamma*T)<v0`.
  Volatility is strictly positive; gamma may be negative. The final reported
  barrier is debt, whereas survival follows the same initial exponential
  barrier formula as R. This is not a separate terminal-debt default jump.
- `at1p(v0,h0,shape,sigma,rate,time)`: `0<h0<v0`, shape B nonnegative;
  sigma/rate vectors have the same length as time. Nonnegative volatilities
  are integrated as `sum(sigma^2*dt)`. The reported firm value uses spot
  `v0*exp(rate*time)` rather than integrating the rate vector, as in R.
- `sbtv(v0,barrier,probability,shape,sigma,rate,time)`: scenario barriers and
  probabilities have equal nonzero length; each barrier lies in `(0,v0)`.
  Nonnegative probabilities must sum to one within `1e-12` (R uses exact
  equality). Survival is their AT1P-weighted mixture; hazards are calculated
  from the mixture, not averaged scenario hazards.
- `cds(time,intensity,rate,premium,recovery,simplified)`: hazard length equals
  dates, hazards are nonnegative, rate length is one or the date count;
  premium is positive and recovery lies in `[0,1)`. With the default false
  flag, piecewise hazards and forward rates drive continuous protection and
  the **source's subtracted accrued-premium** convention. Combined rate zero
  uses analytic integrals. True uses discrete endpoint legs. A zero annuity
  produces an invalid curve; a finite negative annuity is retained to expose
  the source convention, not silently replaced. Price = protection-premium;
  rate = protection/unscaled premium annuity.
- `cds2(time,maturity,rate_time,rate,intensity_time,intensity,premium,recovery,simplified)`:
  interpolate spot rates and hazards linearly onto payment dates, compute
  CDS legs, then interpolate the legs/survival/price onto finite output
  maturities. Output maturities need not be sorted and are clamped outside
  the payment grid. Each input knot grid is finite, sorted and unique;
  singleton grids are supported. The rate is the ratio of interpolated legs.

## Simulation

`merton_sim_from_shocks(v0,rate,sigma,time,shocks)` implements
`V=V0*exp(cumsum((rate-sigma^2/2)*dt+sigma*sqrt(dt)*Z))`. Shocks have shape
`(size(time),number_of_paths)`. It is pure and supports deterministic paired
innovations for R comparison. `merton_sim(v0,rate,sigma,time,paths,seed)` uses
the same formula with a local Park–Miller generator and Box–Muller normal
transformation. It is pure with no global state. Default seed is 104729;
explicit integer seeds are normalized modulo 2147483646 plus one. Its draw
order is column-major and one pair of uniforms is used for each normal,
including at zero-length initial intervals. R's seed, generator, draw
consumption and default wall-clock seed are intentionally not reproduced.

## Calibration

`use creditrisk_calibration` exports four functions, `calibration_result_t`
and the reused dependency's `lbfgsb_control_t`/`lbfgsb_result_t` types.

`calibration_result_t` contains:
- `parameters(:)`: fitted parameters in the order listed below
- `rates(:)`: fitted decimal CDS spreads at quote dates
- `error`: mean squared spread error, independently recomputed at returned parameters
- `initial_error`, `rms_error`, `max_abs_residual`: objective at the starting
  model, final RMS spread residual and largest absolute quote residual
- `warning`: nonempty for unsuccessful optimization or a maximum residual
  above one basis point (0.0001 in decimal spreads); this is a diagnostic
  threshold, not a user-specific financial materiality judgment
- `status`: `credit_ok` when the returned model and objective are valid; this
  does **not** imply optimizer convergence
- `optimizer`: complete L-BFGS-B diagnostics including `success`,
  `convergence`, `message`, iterations and evaluation counts

Check both `status` and `optimizer%success`. Invalid preflight inputs produce
an invalid result; arrays may be unallocated and the objective remains huge.
A valid but iteration-limited fit is returned with its diagnostics for review.
Objective evaluation uses a scaled mean square. Exceptionally large finite
rates/quotes whose magnitudes exceed `sqrt(huge(dp)/400)` are rejected before
subtraction/squaring, leaving representable room for penalty scaling.
There are no global mutable objective data. Optimization is impure because of
the sibling optimizer's implementation.

| Procedure and required arguments | Returned parameter order |
|---|---|
| `calibrate_black_cox(v0,quote,rate,time,...)` | `[L,K,sigma,gamma]` |
| `calibrate_at1p(v0,quote,rate,time,...)` | `[H0,B,sigma(:)]` |
| `calibrate_sbtv(v0,probability,quote,rate,time,...)` | `[B,sigma(:),H(:)]` |
| `calibrate_cds(rate,time,maturity,quote,...)` | hazard values at `maturity` knots |

All allow optional `initial(:)` and `control`. The first, second and fourth
also allow optional `premium,recovery,simplified` before `initial,control`.
BlackCox rate is scalar; AT1P/SBTV rates match time; CDS accepts a singleton
rate vector or payment-grid rates. Market quote vectors match output
maturities. `initial` is in returned parameter order except BlackCox, where
optimizer coordinates are `[L,L-K,sigma,gamma]` to preserve the source
parameterization. The returned optimized L/K fix the R reporting bug.

Default starts exactly follow R: BlackCox `[v0/2,.075*v0,.2,.08]`; AT1P
`[.4,.02,.2,...]`; SBTV `[.02,.2,...,.4*v0,.6*v0,...]`; CDS all `.1`.
AT1P's `.4` is absolute, not scaled by v0. Defaults that violate the physical
barrier domain are rejected; supply a valid initial vector. More than three
SBTV scenarios can make default barriers exceed v0 and likewise require a
custom start. All optimizer coordinates have zero lower bounds. Invalid
trial domains receive finite penalty `max(1,100*initial_error)`, unlike R's possible immediate
nonfinite-objective failure. The objective scale keeps invalid trials expensive
without collapsing line-search steps as a huge constant can do. No arbitrary
data-dependent parameter upper
bounds are added.

Default iteration limits are 100 for BlackCox, 5000 otherwise. Memory=5,
factr=1e7, pgtol=0, finite_difference_step=.001, and custom parameter-change
stopping is disabled. The reused L-BFGS-B step is relative for parameter
magnitudes above one; R uses an absolute .001 step. The optional control
replaces **all** defaults, so construct it deliberately. Optimizer versions,
finite-difference scaling, finite invalid-trial penalties and rounding may
change convergence paths and fitted parameters. These inverse problems may
be weakly identified; no global-optimum guarantee is made.

`calibrate_cds` first interpolates knot hazards to payment dates, then
interpolates computed **rates** to quote dates. Do not substitute `cds2` in
its objective, because ratio-of-interpolated-legs is a different operation.
SBTV's upstream ellipsis is unused: its calibration keeps recovery .4,
premium .005 and the default continuous formula, with no override arguments.

## Supporting helpers

Public `linear_constant(knots,values,targets)` supplies validated linear
interpolation with endpoint clamping, and `valid_times(time,allow_zero)`
validates maturity grids. They are shared by the two package modules and do
not add to the 13-function R coverage count. All maintained Fortran dummy
arguments have explicit intents and trailing FORD documentation.
