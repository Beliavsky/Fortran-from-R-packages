# API quick reference

## High-level portfolio routines

- `parmaspec`: construct a validated portfolio specification.
- `parmasolve`: solve minimum-risk, maximum-reward, risk/reward, or utility
  formulations.
- `parmafrontier`: solve a vector of target returns.
- `parmautility2`, `parmautility4`: CARA moment approximations.
- `spec_risk_value`, `reward_value`, `repair_weights`, `validate_spec`.

## Risk codes

`risk_mad`, `risk_ev`, `risk_minimax`, `risk_cvar`, `risk_cdar`, `risk_lpm`,
`risk_upm`, and `risk_rachev`.

## Objective codes

`solve_min_risk`, `solve_max_reward`, `solve_max_ratio`, and `solve_utility`.

## Standalone solvers

- `cmaes_minimize`
- `qp_box_budget`
- `lp_simplex`
- `milp_binary_solve`
- `socp_solve`

The standalone solver result types are `cmaes_result`, `qp_result`, `lp_result`,
`milp_result`, and `socp_result`.
