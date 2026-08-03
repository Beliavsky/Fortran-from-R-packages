# API reference

All public names are re-exported by module `multi_asset_options`.

## Configuration

### `initialize_config(config, n_asset, status)`

Allocates all per-asset arrays and supplies practical defaults. Set the desired
fields in `config%opt`, `config%fd`, and `config%time` before pricing.

### `validate_config(config, status)`

Checks dimensions, option flags, grid controls, timestepping controls, and the
basic symmetry/range requirements of the correlation matrix.

## Grid and payoff routines

### `node_spacer(strike, left_bound, right_bound, nodes, density, k_shift, x, status)`

Creates the one-dimensional grid used by the upstream `nodeSpacer` routine.

### `build_grid(config, grid, status)`

Creates all asset grids, dimensions, strides, and the total node count.

### `payoff_values(pay_type, pc_flag, strike, grid, values, status)`

Computes the digital, best-of, or worst-of payoff in classical state order.

### `interpolate_value(grid, values, spot, result_value, status)`

Performs multilinear interpolation, clamping points outside the grid to the
nearest boundary.

## Finite-difference operator

### `build_fdm_operator(grid, rf, q, vol, rho, operator, status)`

Constructs the Black-Scholes spatial operator as a native `csr_matrix`.

Useful sparse helpers are `csr_matvec`, `csr_to_dense`, and `csr_diagonal`.

## Pricing

### `price_multi_asset(config, result, status [, solver_tolerance, solver_max_iterations])`

Runs the full backward finite-difference solution. The result contains:

- `grid`: all asset nodes and indexing metadata;
- `value(:,k)`: the option state at time `time(k)`;
- `time`: remaining times to maturity, from maturity to valuation time;
- `linear_iterations`: total sparse iterations at each time step;
- `penalty_iterations`: American penalty iterations at each time step.

## Enumerated integer constants

- Payoffs: `payoff_digital`, `payoff_best_of`, `payoff_worst_of`
- Exercise: `exercise_european`, `exercise_american`
- Direction: `option_call`, `option_put`
- Timesteps: `timestep_constant`, `timestep_adaptive`
