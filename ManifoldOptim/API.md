# API

## Types

### `manifold_component`

Fields correspond to the R manifold definitions: `kind`, `n`, `m`, `p`,
`numofmani`, and `param_set`. Prefer `make_component`.

### `manifold_domain`

Contains `component(:)`. Multiple components form a product manifold.
`domain%length()` returns the flat storage length expected by `x0` and
callbacks.

### `solver_options`

Important fields include:

- `tolerance`, `max_iteration`, `debug`
- `isconvex`
- `memory` (L-BFGS / LRTRSR1 memory)
- `max_linesearch`, `initial_step`, `min_step`, `max_step`
- `line_search`: `LINESEARCH_ARMIJO`, `LINESEARCH_WOLFE`,
  `LINESEARCH_STRONG_WOLFE`, `LINESEARCH_EXACT`, or `LINESEARCH_INPUTFUN`
- `line_search_proc`: procedure pointer used with `LINESEARCH_INPUTFUN`
- `armijo`, `wolfe`
- `qn_nu`, `qn_mu`
- `broyden_phi` (upstream default is 1)
- `eps_numerical_grad`, `eps_numerical_hess`
- `trust_radius`, `max_trust_radius`, `sr1_skip`
- `cg_beta`: `"HS"`, `"FR"`, `"PR"`, `"PR+"`, `"FR-PR"`, `"DY"`, or
  `"HZ"`

### `solver_result`

Contains `xopt`, `fval`, `normgf`, `normgfgf0`, `iter`, evaluation and manifold
operation counts, status/message, and objective/gradient histories.

## Objective and derivative callbacks

```fortran
subroutine objective_callback(x, f)
subroutine gradient_callback(x, g)
subroutine hessvec_callback(x, eta, hess_eta)
```

The analytical gradient and Hessian-vector callbacks are Euclidean derivatives,
as in the upstream ManifoldOptim problem interface. The library converts them
to Riemannian derivatives where the upstream manifold supplies that conversion.

## Optimization

`manifold_optimize` has three forms:

```fortran
call manifold_optimize(domain, x0, objective, method, result [, options])
call manifold_optimize(domain, x0, objective, gradient, method, result [, options])
call manifold_optimize(domain, x0, objective, gradient, hessvec, method, result [, options])
```

Supported method strings are:

```text
RSD RCG RNewton RBFGS LRBFGS RBroydenFamily RWRBFGS
RTRSD RTRNewton RTRSR1 LRTRSR1
```

## Custom line search

The callback signature is:

```fortran
function line_search_callback(x, eta, initial_step, initial_slope) result(step)
  real(dp), intent(in) :: x(:), eta(:), initial_step, initial_slope
  real(dp) :: step
end function
```

Associate it with the option procedure pointer:

```fortran
opt%line_search = LINESEARCH_INPUTFUN
opt%line_search_proc => my_line_search
```

## Manifold utilities

- `project_tangent`
- `euclidean_to_riemannian_gradient`
- `euclidean_hess_to_riemannian`
- `retract_point`
- `differentiated_retraction`
- `transport_vector`
- `inverse_transport_vector`
- `cotangent_vector`
- `manifold_beta`
- `manifold_metric`
- `metric_dual`
- `random_manifold_point`
- `point_is_valid`
- `orthonorm`

## Storage

Matrices are flattened in Fortran column-major order.

LowRank points retain the upstream factor ordering:

```text
vec(U), vec(D), vec(V)
```

with `U` n-by-p, `D` p-by-p, and `V` m-by-p. LowRank tangent vectors use the
same public length, but projection/metric/transport apply the quotient-space
scaling internally.
