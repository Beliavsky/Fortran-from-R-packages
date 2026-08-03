# API

All public facilities are re-exported by:

```fortran
use ldhmm
```

The floating-point kind is `dp = kind(1.0d0)`.

## Main types

### `ecld_type`

Fields:

- `lambda`: lambda shape parameter
- `sigma`: positive scale parameter
- `mu`: location parameter

### `ldhmm_model`

Core fields:

- `m`: number of states
- `param_nbr`: 2 for `(mu,sigma)` or 3 for `(mu,sigma,lambda)`
- `param(m,param_nbr)`: state distribution parameters
- `gamma(m,m)`: row-stochastic transition matrix
- `delta(m)`: initial state probabilities
- `stationary`: whether `delta` is derived from `gamma` during parameter transforms

Fit and decoding results are stored in `mllk`, `aic`, `bic`, `iterations`,
`return_code`, `observations`, `states_prob`, `states_local`, `states_global`,
`states_local_stats`, and `states_global_stats`.

### `ldhmm_fit_control`

Important fields:

- `optimizer`: `bfgs`, `nelder-mead`, or aliases `nlm`, `nm`, `optimx`
- `max_iterations`
- `tolerance`
- `gradient_step`
- `initial_simplex_step`
- `min_gamma`
- `decode`
- `print_level`

## Lambda distribution

```fortran
distribution = ecld_create(lambda, sigma, mu, status)
y = ecld_pdf(distribution, x)
p = ecld_cdf(distribution, x)
q = ecld_ccdf(distribution, x)
r = ecld_random(distribution)
```

`ecld_pdf`, `ecld_cdf`, and `ecld_ccdf` accept either a scalar or a rank-1
array. Analytic statistics are available through `ecld_mean`, `ecld_variance`,
`ecld_sd`, `ecld_skewness`, and `ecld_kurtosis`.

## Model construction and parameters

```fortran
model = ldhmm_create(m, param, gamma_matrix, delta, stationary, optimizer, status)
gamma_matrix = ldhmm_gamma_init(m, p1, p2, probability, min_gamma)
status = ldhmm_validate(model)
working = ldhmm_natural_to_working(model, mu_scale)
model2 = ldhmm_working_to_natural(model, working, mu_scale, status)
call ldhmm_stationary_distribution(gamma_matrix, delta, status)
```

If `delta` is omitted from `ldhmm_create` and `stationary=.true.`, the
stationary distribution is calculated automatically.

## Likelihood and decoding

```fortran
mllk = ldhmm_mllk(model, x, status)
call ldhmm_log_forward(model, x, log_alpha, status)
call ldhmm_log_backward(model, x, log_beta, status)
decoded = ldhmm_decode(model, x, do_global, do_stats, status)
states = ldhmm_viterbi(model, x, status)
```

The forward and backward arrays contain unscaled log probabilities. NaN
observations are treated as missing, with emission likelihood one.

## Conditional distributions and forecasts

```fortran
density = ldhmm_conditional_prob(model, x, grid, status)
state_prob = ldhmm_forecast_state(model, x, horizon, status)
density = ldhmm_forecast_prob(model, x, future_grid, horizon, status)
forecast = ldhmm_forecast_volatility(model, x, future_x, ma_order, days, status)
residuals = ldhmm_pseudo_residuals(model, x, grid_length, status)
```

`ldhmm_forecast_prob` returns `(horizon,size(future_grid))`.
`ldhmm_forecast_state` returns `(m,horizon)`. Volatility forecasts return a
`(2,size(future_x))` array: candidate return and annualized volatility percent.

## Statistics

```fortran
state_stats = ldhmm_ld_stats(model, annualize, days_per_year)
observed_stats = ldhmm_calc_stats_from_obs(model, drop, use_local)
history = ldhmm_decode_stats_history(model, ma_order, annualize, days_per_year)
ma = ldhmm_sma(x, order, na_backfill)
acf = ldhmm_abs_acf(x, lag_max, drop)
clean = ldhmm_drop_outliers(x, drop)
returns = ldhmm_prices_to_log_returns(prices, randomize_zeros, zero_scale)
```

State statistics columns are mean, standard deviation, Pearson kurtosis, and
mean/variance allocation. Observed statistics add skewness, count, and
allocation.

## Simulation

```fortran
simulated = ldhmm_simulate_state_transition(model, init=n, status=status)
advanced = ldhmm_simulate_state_transition(simulated, status=status)
acf = ldhmm_simulate_abs_acf(model, n, lag_max, status)
call seed_random(12345)
```

With `init`, independent initial states are sampled from `delta`. Without
`init`, every state path stored in `states_local` advances one transition.

## Maximum-likelihood fitting

```fortran
control = ldhmm_fit_control()
control%optimizer = 'bfgs'
control%max_iterations = 1000
control%decode = .true.
fitted = ldhmm_fit(model, x, control, status)
```

All positive parameters and transition probabilities are optimized through the
same unconstrained working-parameter representation used by the upstream
package. Gradients for BFGS are central finite differences.

## Status codes

- `LDHMM_SUCCESS`
- `LDHMM_INVALID_ARGUMENT`
- `LDHMM_ALLOCATION_ERROR`
- `LDHMM_NUMERICAL_ERROR`
- `LDHMM_MAX_ITERATIONS`
- `LDHMM_LINE_SEARCH_FAILED`
