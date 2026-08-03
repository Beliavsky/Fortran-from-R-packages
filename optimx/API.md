# API

## Callback

`optimx_callback(x, value, gradient, hessian, need_gradient, need_hessian,
status)` is the single objective callback type. A callback must always return
`value`. It fills `gradient` or `hessian` only when the corresponding flag is
true. Set `status = 0` for a valid evaluation.

## Core types

- `optimx_problem`: callback, dimensions, bounds, and integer mask
- `optimx_control`: iteration/evaluation limits and tolerances
- `optimx_result`: solution, derivatives, counts, status, method, KKT flags
- `optimx_multi_result`: results from several methods or starts
- `derivative_check`, `hessian_check`, `kkt_result`, `bounds_result`

A mask value of zero fixes a parameter. Any nonzero mask value permits it to
move subject to its lower and upper bounds.

## Solver entry points

`optimr` dispatches by method name. Direct entry points are `Rvmmin`,
`Rvmminu`, `Rvmminb`, `Rcgmin`, `Rcgminu`, `Rcgminb`, `nvm`, `ncg`, `hjn`,
`snewton`, `snewtm`, `tn`, and `tnbc`.

`opm` and `optimx` run several methods from the same initial point.
`multistart` treats each column of its start matrix as one initial point.
`polyopt` runs methods sequentially, feeding each result to the next method.
`proptimr` chooses the best converged run.

## Derivatives and diagnostics

- `grfwd`, `grback`, `grcentral`, `grnd`, `grpracma`
- `gHgen`, `gHgenb`
- `fnchk`, `grchk`, `hesschk`
- `kktchk`, `optchk`, `pd_check`, `scalechk`
- `bmchk`, `bmstep`, `axsearch`

## Utilities

`ctrldefault`, `dispdefault`, `checksolver`, `checkallsolvers`, `opm2optimr`,
`optimr2opm`, and the saved context `optsp` are also public.
