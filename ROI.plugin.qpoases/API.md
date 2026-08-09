# API

## Direct solve

```fortran
use qpoases

type(qpoases_result) :: r

call solve_qproblem(H, g, A, lb, ub, lbA, ubA, r, &
                     hessian_type=hst_unknown, max_nwsr=2000)

call solve_qproblemb(H, g, lb, ub, r, &
                      hessian_type=hst_unknown, max_nwsr=2000)
```

`A` is stored in the conventional Fortran shape `(n_constraints,n_variables)`.

`qpoases_result` contains:

* `x` - primal solution.
* `y` - dual solution, ordered as variable bounds followed by general
  constraints.  Positive means a lower-side multiplier and negative an
  upper-side multiplier.
* `objval` - `0.5*x^T*H*x + g^T*x`.
* `status`, `nwsr`.
* state flags `initialised`, `solved`, `infeasible`, `unbounded`.
* working-set counts.

## Persistent model

```fortran
type(qpoases_model) :: model
integer :: status

call init_qproblem(model,H,g,A,lb,ub,lbA,ubA,2000,status)
call hotstart_qproblem(model,g2,lb2,ub2,lbA2,ubA2,2000,status)
```

Simply bounded:

```fortran
call init_qproblemb(model,H,g,lb,ub,2000,status)
call hotstart_qproblemb(model,g2,lb2,ub2,2000,status)
```

SQProblem-style matrix update:

```fortran
call init_sqproblem(model,H,g,A,lb,ub,lbA,ubA,2000,status)
call hotstart_sqproblem(model,H2,g2,A2,lb2,ub2,lbA2,ubA2,2000,status)
```

The model retains its previous primal point and working-set IDs as the warm
start for the next solve.

## Getters

* `get_objval(model)`
* `get_primal_solution(model,x)`
* `get_dual_solution(model,y)`
* `get_number_of_variables(model)`
* `get_number_of_free_variables(model)`
* `get_number_of_fixed_variables(model)`
* `get_number_of_constraints(model)`
* `get_number_of_equality_constraints(model)`
* `get_number_of_active_constraints(model)`
* `get_number_of_inactive_constraints(model)`
* `is_initialised(model)`
* `is_solved(model)`
* `is_infeasible(model)`
* `is_unbounded(model)`

## ROI-style solve

```fortran
use roi_qpoases

character(len=2) :: direction(m)

direction = [">=", "<=", "==", ...]
call roi_solve_qp(Q,L,A,direction,rhs,result)
```

This uses ROI's default variable lower bound of zero if `lower` is omitted.
`maximum=.true.` negates `Q` and `L` for the internal minimizer and reports the
objective in the user's original convention.

## Options

`type(qpoases_options)` mirrors the fields of the R package's
`default_control()`.  Core convergence/regularisation fields are used by the
Fortran solver; advanced qpOASES homotopy-specific fields are retained for API
compatibility but are not all algorithmically active in v0.1.0.
