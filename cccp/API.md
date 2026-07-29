# API reference

All public symbols are re-exported by modules `cccp` and `cccp_api`. Use `cccp_api` to call the original generic function name `cccp`; use `cccp_solve` when importing module `cccp`.

## Kinds and status codes

- `dp`: double-precision real kind.
- `cccp_success`
- `cccp_invalid_input`
- `cccp_infeasible_start`
- `cccp_singular_system`
- `cccp_max_iterations`
- `cccp_domain_error`

## Control and result types

### `type(cccp_control)`

Fields correspond to the R `ctrl()` object:

- `maxiters`
- `max_outer`
- `abstol`
- `reltol`
- `feastol`
- `stepadj`
- `beta`
- `barrier_growth`
- `trace`

`ctrl(...)` constructs this type with optional overrides.

### `type(cccp_solution)`

- `x`: primal variables.
- `y`: equality multipliers.
- `s`: concatenated cone slacks.
- `z`: approximate cone dual variables.
- `cone_offsets`: inclusive ranges for each cone in `s` and `z`.
- `state`: objective, gap, residual, and slack diagnostics.
- `status`: normally `optimal`, `unknown`, `infeasible`, or an error status.
- `niter`: total Newton iterations.
- `info`: integer status code.

Extractors matching the R names are available: `getx`, `gety`, `gets`, `getz`,
`getstate`, `getstatus`, `getniter`, and `getparams`.

## Cone constructors

### `nnoc(G, h)`

Constructs the nonnegative-orthant inequality

```text
G x <= h
```

### `socc(F, g, d, f)`

Constructs the Lorentz-cone inequality

```text
||F x + g||_2 <= d^T x + f
```

The argument ordering and sign convention match the R package.

### `psdc(Flist, F0)`

`Flist(:,:,j)` contains matrix `F_j`. The constraint is

```text
F0 - sum_j x_j F_j  is positive semidefinite.
```

### `nlfc(G, h)`

Compatibility alias of `nnoc`. Nonlinear constraints themselves are supplied by
callback procedures to `dnl` or `dcp`.

## Main solvers

### Linear programs

```fortran
call solve_lp(q, A, b, cones, control, sol)
sol = cccp(q, A, b, cones, control)        ! with use cccp_api
sol = cccp_solve(q, A, b, cones, control)  ! with use cccp
```

Minimizes `q^T x`.

### Quadratic programs

```fortran
call solve_qp(P, q, A, b, cones, control, sol)
sol = cccp(P, q, A, b, cones, control)        ! with use cccp_api
sol = cccp_solve(P, q, A, b, cones, control)  ! with use cccp
```

Minimizes

```text
0.5 x^T P x + q^T x.
```

`P` is symmetrized internally.

### Linear objective with nonlinear inequalities

```fortran
call dnl(q, x0, constraints, mnl, A, b, cones, control, sol)
```

The callback has interface:

```fortran
subroutine constraints(x, f, g, h, info)
   real(dp), intent(in)  :: x(:)
   real(dp), intent(out) :: f(:)
   real(dp), intent(out) :: g(:,:)
   real(dp), intent(out) :: h(:,:,:)
   integer,  intent(out) :: info
end subroutine
```

Here `f(i) <= 0`, `g(i,:)` is its gradient, and `h(:,:,i)` is its Hessian.
The solver performs phase-I recovery when `x0` is in the callback domain but is
not strictly feasible.

### General convex objective

```fortran
call dcp(x0, objective, A, b, cones, control, sol, constraints, mnl)
```

The objective callback has interface:

```fortran
subroutine objective(x, f, g, h, info)
   real(dp), intent(in)  :: x(:)
   real(dp), intent(out) :: f
   real(dp), intent(out) :: g(:)
   real(dp), intent(out) :: h(:,:)
   integer,  intent(out) :: info
end subroutine
```

`solve_dnl` and `solve_dcp` are descriptive aliases of `dnl` and `dcp`.

## Problem-definition compatibility types

- `dlp(q, A, b, cones)` returns `type(dlp_problem)`.
- `dqp(P, q, A, b, cones)` returns `type(dqp_problem)`.
- `cps(problem, control)` solves either type.

This preserves the original define-then-solve workflow without R reference
classes.

## Specialized algorithms

### `l1(P, q, control)`

Solves

```text
minimize ||P u - q||_1
```

through the same linear-cone solver. The first `size(P,2)` elements of `x` are
`u`; the remaining elements are absolute-residual epigraph variables.

### `rp(x0, P, mrc, control)`

Computes a long-only risk-parity portfolio by minimizing the package's convex
log-barrier objective and normalizing the solution to sum to one. `mrc` is
normalized automatically.

### `gp(F0, g0, constraints, nno, A, b, control)`

Solves a geometric program in logarithmic variables. A posynomial is represented
by `type(gp_function)` with fields:

- `f`: exponent matrix.
- `g`: logarithms of coefficients.

The returned `x` is transformed back to the positive original scale.

### `logsumexp_value_gradient_hessian`

Evaluates the log-sum-exp value, analytical gradient, and analytical Hessian used
by geometric programming.

## Cone algebra

The following low-level routines expose the computational cone operations used
for testing and specialized numerical work:

- `cone_inner_product`
- `cone_norm`
- `cone_identity`
- `cone_jordan_product`
- `cone_jordan_inverse`
- `cone_max_step`
- `cones_interior`
- `minimum_cone_slack`

They support nonnegative, Lorentz, and positive-semidefinite cones.
