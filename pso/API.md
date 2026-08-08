# API

The main public module is:

```fortran
use pso
```

## `psoptim`

Vector bounds:

```fortran
call psoptim(par, fn, lower, upper, result [, control] [, gr])
```

Scalar bounds:

```fortran
call psoptim(par, fn, lower_scalar, upper_scalar, result [, control] [, gr])
```

`par` determines problem dimension.  If every entry is finite/non-NaN and lies
inside the bounds, the first particle starts at `par`.  Otherwise use IEEE NaNs
to request fully random initialization, matching the role of `NA` in R.

Objective callback:

```fortran
function fn(x) result(f)
   use pso, only : dp
   real(dp), intent(in) :: x(:)
   real(dp) :: f
end function fn
```

Optional gradient callback used by hybrid local refinement:

```fortran
subroutine gr(x, g)
   use pso, only : dp
   real(dp), intent(in) :: x(:)
   real(dp), intent(out) :: g(:)
end subroutine gr
```

## `pso_control`

Important fields and R counterparts:

| Fortran | R `control` | Default |
|---|---|---:|
| `trace` | `trace > 0` | `.false.` |
| `fnscale` | `fnscale` | `1` |
| `maxit` | `maxit` | `1000` |
| `maxf` | `maxf` | effectively unlimited |
| `abstol` | `abstol` | `-huge()` |
| `reltol` | `reltol` | `0` |
| `report` | `REPORT` | `10` |
| `swarm_size` | `s` | automatic when `0` |
| `k` | `k` | `3` |
| `informant_p` | `p` | automatic when `<0` |
| `w0`, `w1` | `w` | `1/(2*log(2))` |
| `c_p` | `c.p` | `0.5+log(2)` |
| `c_g` | `c.g` | `0.5+log(2)` |
| `diameter` | `d` | automatic when `<0` |
| `v_max` | `v.max` | disabled when `<0` |
| `rand_order` | `rand.order` | `.true.` |
| `max_restart` | `max.restart` | effectively unlimited |
| `maxit_stagnate` | `maxit.stagnate` | effectively unlimited |
| `vectorize` | `vectorize` | `.false.` |
| `hybrid` | `hybrid` | `pso_hybrid_off` |
| `trace_stats` | `trace.stats` | `.false.` |
| `pso_type` | `type` | `pso_spso2007` |

PSO constants:

```fortran
pso_spso2007
pso_spso2011
```

Hybrid constants:

```fortran
pso_hybrid_off
pso_hybrid_on
pso_hybrid_improved
```

The self-contained local solver additionally exposes these controls through
`pso_control`:

```fortran
hybrid_maxit
hybrid_memory
hybrid_reltol
```

## `pso_result`

Main fields:

```text
par
value
function_evaluations
iterations
restarts
convergence
message
```

Convergence codes follow upstream:

- `0`: absolute tolerance reached
- `1`: function-evaluation limit reached
- `2`: iteration limit reached
- `3`: restart limit reached
- `4`: stagnation limit reached

When both `trace` and `trace_stats` are enabled, these are populated:

```text
ntrace
trace_it
trace_error
trace_f
trace_x
```

`trace_f(:,k)` contains particle fitness values at report `k`, while
`trace_x(:,:,k)` contains the corresponding swarm positions.

## Random seeding

```fortran
call seed_random(integer_seed)
```

This initializes the compiler's intrinsic `random_number` generator
reproducibly for that compiler/runtime.  It is not intended to reproduce R's
Mersenne-Twister stream bit-for-bit.

## Benchmark support

Construct one of the upstream benchmark problems:

```fortran
call make_test_problem("rastrigin", problem)
```

Available names:

```text
parabola
griewank
rosenbrock
rastrigin
ackley
```

Run repeated tests:

```fortran
call run_test_problem(problem, summary [, control])
```

Calculate the empirical success curve:

```fortran
call get_success_curve(summary, feval, rate)
```

Calculate the upstream efficiency statistic:

```fortran
eff = test_efficiency(summary)
```

or collect objective mean/sample standard deviation/minimum/maximum, success
rate, efficiency, and timing:

```fortran
stats = summarize_test(summary)
```
