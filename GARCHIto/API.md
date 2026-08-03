# API

All public entities are provided by module `garchito`.

## Kinds and status constants

```fortran
integer, parameter :: dp
integer, parameter :: garchito_success = 0
integer, parameter :: garchito_max_iterations = 1
integer, parameter :: garchito_invalid_input = 2
integer, parameter :: garchito_numerical_failure = 3
```

## Control type

```fortran
type :: garchito_control
   integer :: max_iterations = 5000
   integer :: max_evaluations = 100000
   real(dp) :: tolerance = 1.0e-9_dp
   real(dp) :: simplex_scale = 0.10_dp
   integer :: trace = 0
end type
```

`trace=0` suppresses optimizer output. A positive value prints the current
objective every `trace` iterations.

## Result type

```fortran
type :: garchito_result
   real(dp), allocatable :: coefficients(:)
   character(len=16), allocatable :: coefficient_names(:)
   real(dp), allocatable :: sigma(:)
   real(dp) :: pred
   real(dp) :: objective
   integer :: convergence
   integer :: iterations
   integer :: evaluations
   character(len=160) :: message
end type
```

`sigma` contains fitted conditional integrated variances, matching the R
package's field name. `pred` is the next conditional integrated variance.

## Unified GARCH-Ito

```fortran
subroutine unified_est(rv, returns, result, control)
   real(dp), intent(in) :: rv(:), returns(:)
   type(garchito_result), intent(out) :: result
   type(garchito_control), intent(in), optional :: control
end subroutine
```

The coefficient order is `omega_g`, `beta_g`, `gamma`. The recursion is

```text
h(t) = omega_g + gamma*h(t-1) + beta_g*returns(t-1)^2
```

with `h(1)=omega_g/(1-beta_g-gamma)`.

## Realized GARCH-Ito

```fortran
subroutine realized_est(rv, result, jv, control)
   real(dp), intent(in) :: rv(:)
   type(garchito_result), intent(out) :: result
   real(dp), intent(in), optional :: jv(:)
   type(garchito_control), intent(in), optional :: control
end subroutine
```

Without `jv`, the coefficient order is `omega_g`, `alpha_g`, `gamma` and

```text
h(t) = omega_g + gamma*h(t-1) + alpha_g*rv(t-1)
```

With `jv`, the coefficient order is `omega_g`, `alpha_g`, `beta_g`, `gamma`
and

```text
h(t) = omega_g + gamma*h(t-1) + alpha_g*rv(t-1) + beta_g*jv(t-1)
```

The jump model uses the sample median of `jv` in its unconditional initial
variance, as in the R source.

## Realized GARCH-Ito with option volatility

```fortran
subroutine realized_est_option(rv, nv, result, jv, homogeneous, control)
   real(dp), intent(in) :: rv(:), nv(:)
   type(garchito_result), intent(out) :: result
   real(dp), intent(in), optional :: jv(:)
   logical, intent(in), optional :: homogeneous
   type(garchito_control), intent(in), optional :: control
end subroutine
```

`homogeneous` defaults to `.true.`. The measurement equation is

```text
nv(t) = b + a*h(t) + error(t)
```

The homogeneous model uses constant error variance `sigma_e^2`. The
heterogeneous model uses `sigma_e^2*h(t)^zeta`.

Coefficient orders:

- no jumps, homogeneous: `omega_g, alpha_g, gamma, a, b, sigma_e`
- no jumps, heterogeneous: the above plus `zeta`
- jumps, homogeneous: `omega_g, alpha_g, beta_g, gamma, a, b, sigma_e`
- jumps, heterogeneous: the above plus `zeta`

## Input rules

- All arrays must be finite.
- `rv` and optional `jv` must be nonnegative.
- Paired series must have equal lengths.
- At least two observations are required for non-option models and at least
  three for the option model.
- A best finite estimate is returned when the optimizer reaches its iteration
  limit; `convergence` and `message` identify that case.
