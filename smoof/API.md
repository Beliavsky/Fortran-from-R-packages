# API summary

All real arithmetic uses `dp = kind(1.0d0)` from `smoof_kinds`.

## Single objective

Representative calls:

```fortran
use smoof_single, only : ackley, hartmann, shekel, michalewicz
real(dp) :: value
value = ackley(x)
value = hartmann(x)             ! length 3, 4, or 6
value = shekel(x4, 10)
value = michalewicz(x, 10.0_dp)
```

Parameterized functions use optional or explicit scalar/vector arguments where
R used generator closures, e.g. `deflected_corrugated_spring(x,k,alpha)` and
`modified_rastrigin(x,k)`.

`swiler2014(x1,x2,x3)` uses integer `x1=1..5` plus two real arguments.
The upstream Inverted Vincent benchmark is `vincent(x)`; its metadata indicates
maximization, which is an optimizer concern rather than a different formula.

## Multi-objective

```fortran
call dtlz2(x, n_objectives, f)
call wfg1(z, n_objectives, k, f)
call zdt1(x, f2)
call mop5(x, f3)
call uf10(x, f3)
call mmf15a(x, n_objectives, np, f)
```

## NK landscapes

`nk_landscape` stores flat arrays rather than R lists. `evaluate()` reproduces
the upstream native offset calculation and returns the mean contribution.
`nk_evaluate_raw` is available when callers already own flat arrays.
