# multiAssetOptions-fortran

A self-contained modern Fortran/FPM translation of the computational core of
R package `multiAssetOptions` 0.1-2.

The library values European and American options on one or more assets under a
constant-volatility, constant-correlation Black-Scholes model. It supports
nonuniform spatial grids, mixed derivative terms, theta timestepping,
Rannacher smoothing, penalty iteration, and adaptive timesteps.

## Implemented functionality

- Uniform and strike-concentrated node spacing
- Optional mesh shifting to place a strike on or between nodes
- Digital, best-of, and worst-of payoffs with call/put flags per asset
- Native CSR finite-difference operator construction
- Nonuniform first and second derivatives
- Correlated mixed derivatives with inward first-order boundary treatment
- European and American exercise
- Theta timestepping and Rannacher smoothing
- Forsyth-Vetzal-style adaptive timesteps
- Penalty iteration for the American exercise constraint
- Jacobi-preconditioned BiCGSTAB sparse linear solves
- Classical state ordering with the first asset varying fastest
- Multilinear interpolation at a supplied spot vector

Plotting and R-specific list, S3/S4, data-frame, and animation infrastructure
are intentionally omitted.

## Build

```sh
fpm build
fpm test
fpm run --example european_call
fpm run demo_multi_asset_options
```

No external numerical library is required.

## Minimal example

```fortran
program price_call
   use multi_asset_options
   implicit none

   type(pricing_config) :: config
   type(pricing_result) :: result
   type(status_type) :: status
   real(dp) :: value

   call initialize_config(config,1,status)
   if (.not. status%ok()) error stop status%message

   config%opt%pay_type = payoff_best_of
   config%opt%exercise_type = exercise_european
   config%opt%pc_flag = [option_call]
   config%opt%strike = [100.0_dp]
   config%opt%ttm = 1.0_dp
   config%opt%rf = 0.05_dp
   config%opt%q = [0.0_dp]
   config%opt%vol = [0.2_dp]
   config%fd%m = [120]
   config%fd%k_mult = [4.0_dp]
   config%fd%k_shift = [2]
   config%time%n_steps = 240

   call price_multi_asset(config,result,status)
   if (.not. status%ok()) error stop status%message

   call interpolate_value(result%grid, &
      result%value(:,size(result%time)),[100.0_dp],value,status)
   if (.not. status%ok()) error stop status%message

   print *, value
end program price_call
```

## Result layout

`result%value(:,1)` is the payoff at maturity. Later columns move backward in
calendar time until `result%time(end) = 0`. For multiple assets, the first
asset index varies fastest, matching the classical ordering documented by the
R package.

The number of spatial nodes grows as the product of all per-asset node counts.
High-dimensional problems therefore remain subject to the usual finite-
difference curse of dimensionality.

## License

`GPL-2.0-only OR GPL-3.0-only`. See `LICENSE` and `licenses/`.
