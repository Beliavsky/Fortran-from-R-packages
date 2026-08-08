# Porting notes

## Data model

R lists, closures, environments, and dynamically assembled stage/sub-stage
objects are replaced by `mize_control_t`, `mize_state_t`, and `mize_result_t`.
The state object is public so an application can inspect or persist an
optimizer between calls to `mize_step`.

## Callbacks

The central callback returns both objective and gradient. Optional callbacks
provide a Hessian, Hessian-vector product, monitor/cancellation hook, custom
momentum schedule, and polymorphic user data. Every procedure dummy has an
explicit abstract interface; the release is compiled with
`-Werror=implicit-interface`.

## Line searches

The R package contains separate long implementations of More-Thuente,
Rasmussen, Schmidt/minFunc, and Hager-Zhang searches. This port preserves their
public names and configurable Armijo/weak/strong-Wolfe conditions but routes
those names through one native safeguarded interpolation-and-zoom Wolfe engine.
Consequently accepted iterates and evaluation counts need not match R exactly.
Constant, backtracking, and bold-driver methods are separate implementations.

## Newton information

Newton and partial-Hessian methods accept a dense Hessian callback. Without
one, a symmetric finite-difference Hessian is formed from gradient evaluations.
Truncated Newton accepts a Hessian-vector callback, a dense Hessian callback,
or a finite-difference product.

## Momentum

Built-in constant, switch, ramp, and Nesterov-convex schedules are available.
A custom schedule can be supplied as a procedure argument. Function- and
gradient-based restart reset momentum for subsequent iterations.

## Numerical choices

- BFGS/L-BFGS updates are skipped when curvature is too small.
- SR1 uses the standard denominator safeguard and resets to identity if its
  approximation ceases to define a descent direction.
- Newton directions fall back to steepest descent if the dense solve fails or
  does not produce descent.
- All trial objective and gradient values are checked for IEEE finiteness.
- Dense linear solves use self-contained partial-pivoting elimination; no BLAS
  or LAPACK dependency is required.

## Omitted R infrastructure

R documentation generation, vignettes, list formatting, S3-like printing,
lifecycle/aspect hooks, and plotting are not translated. The original source
is bundled for provenance.
