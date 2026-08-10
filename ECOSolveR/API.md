# API

Use the umbrella module:

```fortran
use ecos
```

## Core problem

The mathematical problem is

```text
minimize    c^T x
subject to  A x = b
            h - G x in K
```

with `K` consisting of a nonnegative orthant, zero or more second-order cones, and zero or more exponential cones.

## `type(ecos_dims)`

```fortran
integer :: l = 0
integer, allocatable :: q(:)
integer :: e = 0
```

Rows of `G`/`h` are ordered as `l` positive-orthant rows, each SOC block in `q`, then `e` exponential-cone triples.

## `type(ecos_settings)`

Important sparse controls in v0.4:

```fortran
logical  :: sparse_kkt = .true.
logical  :: sparse_amd = .true.
logical  :: sparse_rcm = .false.

logical  :: equilibrate = .true.
integer  :: equilibration_iters = 3
real(dp) :: equilibration_min = 1.0e-4_dp
real(dp) :: equilibration_max = 1.0e4_dp

logical  :: dynamic_regularization = .true.
integer  :: max_regularization_updates = 6
real(dp) :: regularization = 1.0e-9_dp

integer  :: iterative_refinement = 3
real(dp) :: refinement_tol = 1.0e-10_dp

integer  :: certificate_maxit = 80
integer  :: dense_diagnostic_limit = 256
```

`ecos_control(...)` accepts the commonly adjusted subset, including:

```fortran
ctrl = ecos_control( &
    maxit=100, &
    sparse_kkt=.true., &
    sparse_amd=.true., &
    sparse_rcm=.false., &
    equilibrate=.true., &
    equilibration_iters=3, &
    dynamic_regularization=.true., &
    max_regularization_updates=6, &
    iterative_refinement=3, &
    refinement_tol=1.0e-10_dp, &
    certificate_maxit=80)
```

Set both `sparse_amd=.false.` and `sparse_rcm=.false.` to use identity ordering. If both AMD and RCM are true, AMD takes precedence.

## `type(ecos_result)`

The standard fields include:

```fortran
real(dp), allocatable :: x(:), y(:), s(:), z(:)
integer :: exitflag, iter, mi_iter
real(dp) :: pcost, dcost, pres, dres, gap, relgap
character(len=96) :: infostring
```

Sparse diagnostics include:

```fortran
logical :: sparse_backend_used
integer :: kkt_nnz
integer :: ldl_nnz
integer :: iterative_refinements
integer :: symbolic_analyses
integer :: numeric_factorizations
integer :: regularization_updates
integer :: cached_symbolic_reuses
integer :: cached_warm_starts
integer :: bb_symbolic_reuses
real(dp) :: factor_fill_ratio
real(dp) :: time_ordering
real(dp) :: time_factorization
real(dp) :: time_refinement
real(dp) :: min_col_scale, max_col_scale
real(dp) :: min_row_scale, max_row_scale
```

Sparse certificate output:

```fortran
logical :: primal_certificate_valid
logical :: dual_certificate_valid
real(dp), allocatable :: primal_certificate(:)
real(dp), allocatable :: dual_certificate(:)
```

For `ECOS_PINF`, `primal_certificate` contains the normalized dual ray `(y,z)`. For `ECOS_DINF`, `dual_certificate` contains the normalized primal ray `d`.

## Sparse matrices

`type(ecos_csc_matrix)` uses one-based CSC storage:

```fortran
integer :: nrow, ncol
integer, allocatable :: colptr(:)
integer, allocatable :: rowind(:)
real(dp), allocatable :: values(:)
```

`type(ecos_csr_matrix)` is the row-oriented companion type.

Helpers:

```fortran
call make_csc_matrix(dense,csc)
call csc_from_zero_based(nrow,ncol,matbeg,matind,values,csc,ierr)
```

## One-shot solve

Dense:

```fortran
call ecos_csolve(c,G,h,dims,result,A,b, &
                 bool_vars,int_vars,control,ierr)
```

Sparse:

```fortran
call ecos_csolve(c,Gcsc,h,dims,result,Acsc,b,control=ctrl)
```

For a successful CSC continuous solve with `sparse_kkt=.true.`, the main Newton path does not convert `G` or `A` to dense form.

## Persistent workspace and reuse

```fortran
call setup_problem_csc(problem,c,Gcsc,h,dims,Acsc,b)
call ecos_setup(workspace,problem,control)

call ecos_solve(workspace,result1)

call ecos_update(workspace,c=new_c,h=new_h,b=new_b)
call ecos_solve(workspace,result2)
```

The second solve can reuse both symbolic KKT analysis and the previous interior point.

A same-pattern sparse matrix value update also retains symbolic analysis:

```fortran
call ecos_update(workspace,g_csc=Gnew)
call ecos_solve(workspace,result3)
```

If `Gnew` has the same `colptr`/`rowind` pattern, `result3%cached_symbolic_reuses` is incremented and no new symbolic analysis is required. The previous warm point is invalidated because matrix values and equilibration can change.

A structural pattern change invalidates both appropriately.

## Inaccurate status

The status constants include:

```text
ECOS_OPTIMAL       =  0
ECOS_PINF          =  1
ECOS_DINF          =  2
ECOS_INACC_OFFSET  = 10
ECOS_MAXIT         = -1
ECOS_NUMERICS      = -2
ECOS_OUTCONE       = -3
ECOS_SIGINT        = -4
ECOS_FATAL         = -7
```

If the best restored sparse iterate misses exact tolerances but satisfies `feastol_inacc` and the inaccurate gap tolerances, the solver returns:

```text
ECOS_OPTIMAL + ECOS_INACC_OFFSET
```

## Sparse ECOS_BB

Integer/boolean variables are supplied exactly as before:

```fortran
call ecos_csolve(c,Gcsc,h,dims,result, &
                 bool_vars=bool_index,int_vars=int_index)
```

Sparse node problems retain sparse storage. `result%bb_symbolic_reuses` reports symbolic factorization reuse across branch-and-bound nodes.

## MatrixExtra adapter

The optional project under `integration/matrixextra-adapter/` provides:

```fortran
ecos_csc_from_matrix
ecos_csc_from_coo
setup_problem_matrixextra
```

It copies sparse storage arrays without an intermediate dense matrix.

## Exponential cone convention

For each slack triple `(a,b,c)`:

```text
b > 0, c > 0, a <= c*log(b/c)
```

which is equivalent to the ECOS convention `exp(a/c) <= b/c`, `c > 0`.
