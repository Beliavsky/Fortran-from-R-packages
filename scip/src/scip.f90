! Apache-2.0
!
! Modern Fortran interface corresponding to the computational API of the
! R package scip 1.10.0-3.  Optimization is performed by the package's
! vendored SCIP 10.0.2 / SoPlex 8.0.2 backend through a small plain-C ABI.
module scip
    use, intrinsic :: iso_c_binding
    use, intrinsic :: ieee_arithmetic
    use scip_c_api
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)
    integer, parameter :: PARAM_LOGICAL = 1, PARAM_INTEGER = 2
    integer, parameter :: PARAM_REAL = 3, PARAM_TEXT = 4

    type, public :: scip_info
        real(dp) :: solve_time = 0.0_dp
        integer(c_long_long) :: nodes = 0_c_long_long
        integer(c_long_long) :: iterations = 0_c_long_long
        real(dp) :: gap = 0.0_dp
        integer :: sol_count = 0
    end type

    type, public :: scip_solution
        real(dp) :: objval = 0.0_dp
        real(dp), allocatable :: x(:)
        logical :: available = .false.
    end type

    type, public :: scip_result
        character(len=32) :: status = "unknown"
        real(dp) :: objval = 0.0_dp
        real(dp), allocatable :: x(:)
        integer :: sol_count = 0
        real(dp) :: gap = 0.0_dp
        type(scip_info) :: info
    end type

    type, public :: scip_csc_matrix
        integer :: nrow = 0
        integer :: ncol = 0
        integer, allocatable :: row(:)
        integer, allocatable :: colptr(:)
        real(dp), allocatable :: val(:)
    end type

    type :: scip_parameter
        character(len=:), allocatable :: name
        integer :: kind = 0
        logical :: lvalue = .false.
        integer(c_long_long) :: ivalue = 0_c_long_long
        real(dp) :: rvalue = 0.0_dp
        character(len=:), allocatable :: tvalue
    end type

    type, public :: scip_control
        logical :: verbose = .true.
        integer :: verbosity_level = 3
        integer :: display_freq = 100
        real(dp) :: time_limit = huge(1.0_dp)
        integer(c_long_long) :: node_limit = -1_c_long_long
        integer(c_long_long) :: stall_node_limit = -1_c_long_long
        integer :: sol_limit = -1
        integer :: best_sol_limit = -1
        real(dp) :: mem_limit = huge(1.0_dp)
        integer :: restart_limit = -1
        real(dp) :: gap_limit = 0.0_dp
        real(dp) :: abs_gap_limit = 0.0_dp
        real(dp) :: feastol = 1.0e-6_dp
        real(dp) :: dualfeastol = 1.0e-7_dp
        real(dp) :: epsilon = 1.0e-9_dp
        logical :: presolving = .true.
        integer :: presolve_rounds = -1
        integer :: lp_threads = 1
        integer(c_long_long) :: lp_iteration_limit = -1_c_long_long
        logical :: lp_scaling = .true.
        character(len=1) :: branching_score = 'p'
        character(len=10) :: heuristics_emphasis = 'default'
        integer :: threads = 1
        type(scip_parameter), allocatable, private :: extra(:)
    contains
        procedure :: add_param_logical => control_add_param_logical
        procedure :: add_param_integer => control_add_param_integer
        procedure :: add_param_int64 => control_add_param_int64
        procedure :: add_param_real => control_add_param_real
        procedure :: add_param_text => control_add_param_text
        generic :: add_param => add_param_logical, add_param_integer, add_param_int64, &
                                add_param_real, add_param_text
    end type

    type, public :: scip_model_t
        private
        type(c_ptr) :: handle = c_null_ptr
    contains
        procedure :: add_var => model_add_var
        procedure :: add_vars => model_add_vars
        procedure :: add_linear_cons => model_add_linear_cons
        procedure :: add_quadratic_cons => model_add_quadratic_cons
        procedure :: add_sos1_cons => model_add_sos1_cons
        procedure :: add_sos2_cons => model_add_sos2_cons
        procedure :: add_indicator_cons => model_add_indicator_cons
        procedure :: set_param_integer => model_set_param_integer
        procedure :: set_param_int64 => model_set_param_int64
        procedure :: set_param_real => model_set_param_real
        procedure :: set_param_logical => model_set_param_logical
        procedure :: set_param_text => model_set_param_text
        generic :: set_param => set_param_integer, set_param_int64, &
                                set_param_real, set_param_logical, set_param_text
        procedure :: set_objective_sense => model_set_objective_sense
        procedure :: optimize => model_optimize
        procedure :: get_status => model_get_status
        procedure :: get_solution => model_get_solution
        procedure :: get_objval => model_get_objval
        procedure :: get_nsols => model_get_nsols
        procedure :: get_sol => model_get_sol
        procedure :: get_info => model_get_info
        procedure :: get_nvars => model_get_nvars
        procedure :: get_nconss => model_get_nconss
        procedure :: free => model_free
    end type

    public :: scip_model, scip_solve, make_csc_matrix
    public :: scip_add_var, scip_add_vars, scip_add_linear_cons
    public :: scip_add_quadratic_cons, scip_add_sos1_cons, scip_add_sos2_cons
    public :: scip_add_indicator_cons, scip_set_param, scip_set_objective_sense
    public :: scip_optimize, scip_get_status, scip_get_solution, scip_get_objval
    public :: scip_get_nsols, scip_get_sol, scip_get_info, scip_model_free

    interface scip_solve
        module procedure scip_solve_dense
        module procedure scip_solve_csc
    end interface
    interface scip_add_var
        module procedure model_add_var
    end interface
    interface scip_add_vars
        module procedure model_add_vars
    end interface
    interface scip_add_linear_cons
        module procedure model_add_linear_cons
    end interface
    interface scip_add_quadratic_cons
        module procedure model_add_quadratic_cons
    end interface
    interface scip_add_sos1_cons
        module procedure model_add_sos1_cons
    end interface
    interface scip_add_sos2_cons
        module procedure model_add_sos2_cons
    end interface
    interface scip_add_indicator_cons
        module procedure model_add_indicator_cons
    end interface
    interface scip_set_param
        module procedure model_set_param_integer
        module procedure model_set_param_int64
        module procedure model_set_param_real
        module procedure model_set_param_logical
        module procedure model_set_param_text
    end interface
    interface scip_set_objective_sense
        module procedure model_set_objective_sense
    end interface
    interface scip_optimize
        module procedure model_optimize
    end interface
    interface scip_get_status
        module procedure model_get_status
    end interface
    interface scip_get_solution
        module procedure model_get_solution
    end interface
    interface scip_get_objval
        module procedure model_get_objval
    end interface
    interface scip_get_nsols
        module procedure model_get_nsols
    end interface
    interface scip_get_sol
        module procedure model_get_sol
    end interface
    interface scip_get_info
        module procedure model_get_info
    end interface
    interface scip_model_free
        module procedure model_free
    end interface

contains

    function c_string(s) result(cs)
        character(len=*), intent(in) :: s
        character(c_char), allocatable :: cs(:)
        integer :: i, n
        n = len_trim(s)
        allocate(cs(n + 1))
        do i = 1, n
            cs(i) = s(i:i)
        end do
        cs(n + 1) = c_null_char
    end function

    subroutine make_c_string(s, cs)
        character(len=*), intent(in) :: s
        character(c_char), allocatable, intent(out) :: cs(:)
        integer :: i, n
        n = len_trim(s)
        allocate(cs(n + 1))
        do i = 1, n
            cs(i) = s(i:i)
        end do
        cs(n + 1) = c_null_char
    end subroutine

    function f_string(cs) result(s)
        character(c_char), intent(in) :: cs(:)
        character(len=32) :: s
        integer :: i, n
        s = ''
        n = min(size(cs), len(s))
        do i = 1, n
            if (cs(i) == c_null_char) exit
            s(i:i) = cs(i)
        end do
    end function

    subroutine return_error(rc, where, ierr)
        integer(c_int), intent(in) :: rc
        character(len=*), intent(in) :: where
        integer, intent(out), optional :: ierr
        character(len=160) :: msg
        if (present(ierr)) then
            ierr = int(rc)
        else if (rc /= 0_c_int) then
            write(msg, '(a,": backend error ",i0)') trim(where), int(rc)
            error stop trim(msg)
        end if
    end subroutine

    subroutine validate_model(self)
        class(scip_model_t), intent(in) :: self
        if (.not. c_associated(self%handle)) error stop 'SCIP model has been freed'
    end subroutine

    function scip_model(name, ierr) result(model)
        character(len=*), intent(in), optional :: name
        integer, intent(out), optional :: ierr
        type(scip_model_t) :: model
        character(c_char), allocatable :: cname(:)
        integer(c_int) :: rc
        if (present(name)) then
            call make_c_string(name, cname)
        else
            call make_c_string('scip_model', cname)
        end if
        rc = c_scipf_model_create(cname, model%handle)
        call return_error(rc, 'scip_model', ierr)
    end function

    subroutine model_free(self, ierr)
        class(scip_model_t), intent(inout) :: self
        integer, intent(out), optional :: ierr
        integer(c_int) :: rc
        if (.not. c_associated(self%handle)) then
            if (present(ierr)) ierr = 0
            return
        end if
        rc = c_scipf_model_free(self%handle)
        self%handle = c_null_ptr
        call return_error(rc, 'scip_model_free', ierr)
    end subroutine

    integer function model_add_var(self, obj, lb, ub, vtype, name, ierr) result(index1)
        class(scip_model_t), intent(inout) :: self
        real(dp), intent(in) :: obj
        real(dp), intent(in), optional :: lb, ub
        character(len=*), intent(in), optional :: vtype, name
        integer, intent(out), optional :: ierr
        real(dp) :: l, u
        character(len=1) :: vt
        character(len=64) :: nm
        character(c_char), allocatable :: cnm(:)
        integer(c_int) :: rc, idx
        integer :: nv
        call validate_model(self)
        l = 0.0_dp
        u = ieee_value(0.0_dp, ieee_positive_inf)
        if (present(lb)) l = lb
        if (present(ub)) u = ub
        vt = 'C'
        if (present(vtype)) vt = vtype(1:1)
        if (vt == 'B' .or. vt == 'b') then
            l = max(l, 0.0_dp)
            u = min(u, 1.0_dp)
        end if
        if (present(name)) then
            nm = name
        else
            nv = self%get_nvars()
            write(nm, '("x",i0)') nv + 1
        end if
        cnm = c_string(trim(nm))
        rc = c_scipf_add_var(self%handle, obj, l, u, int(iachar(vt), c_int), cnm, idx)
        call return_error(rc, 'scip_add_var', ierr)
        if (rc == 0_c_int) then
            index1 = int(idx)
        else
            index1 = -1
        end if
    end function

    integer function model_add_vars(self, obj, lb, ub, vtype, names, ierr) result(first_index)
        class(scip_model_t), intent(inout) :: self
        real(dp), intent(in) :: obj(:)
        real(dp), intent(in), optional :: lb(:), ub(:)
        character(len=*), intent(in), optional :: vtype(:), names(:)
        integer, intent(out), optional :: ierr
        integer :: j, idx, er, n0
        real(dp) :: l, u
        character(len=1) :: vt
        character(len=64) :: nm
        if (present(lb)) then
            if (size(lb) /= size(obj)) error stop 'scip_add_vars: lb size mismatch'
        end if
        if (present(ub)) then
            if (size(ub) /= size(obj)) error stop 'scip_add_vars: ub size mismatch'
        end if
        if (present(vtype)) then
            if (size(vtype) /= size(obj)) error stop 'scip_add_vars: vtype size mismatch'
        end if
        if (present(names)) then
            if (size(names) /= size(obj)) error stop 'scip_add_vars: names size mismatch'
        end if
        first_index = self%get_nvars() + 1
        n0 = self%get_nvars()
        do j = 1, size(obj)
            l = 0.0_dp
            u = ieee_value(0.0_dp, ieee_positive_inf)
            vt = 'C'
            if (present(lb)) l = lb(j)
            if (present(ub)) u = ub(j)
            if (present(vtype)) vt = vtype(j)(1:1)
            if (present(names)) then
                idx = self%add_var(obj(j), l, u, vt, names(j), er)
            else
                write(nm, '("x",i0)') n0 + j
                idx = self%add_var(obj(j), l, u, vt, trim(nm), er)
            end if
            if (er /= 0) then
                if (present(ierr)) ierr = er
                first_index = -1
                return
            end if
        end do
        if (present(ierr)) ierr = 0
    end function

    integer function model_add_linear_cons(self, vars, coefs, lhs, rhs, name, ierr) &
        result(index1)
        class(scip_model_t), intent(inout) :: self
        integer, intent(in) :: vars(:)
        real(dp), intent(in) :: coefs(:)
        real(dp), intent(in), optional :: lhs, rhs
        character(len=*), intent(in), optional :: name
        integer, intent(out), optional :: ierr
        integer(c_int), allocatable :: cv(:)
        real(c_double), allocatable :: cc(:)
        real(dp) :: l, r
        character(len=64) :: nm
        character(c_char), allocatable :: cnm(:)
        integer(c_int) :: rc, idx
        integer :: n, nc
        call validate_model(self)
        n = size(vars)
        if (size(coefs) /= n) error stop 'scip_add_linear_cons: size mismatch'
        allocate(cv(max(1, n)), cc(max(1, n)))
        if (n > 0) then
            cv(1:n) = int(vars, c_int)
            cc(1:n) = real(coefs, c_double)
        end if
        l = ieee_value(0.0_dp, ieee_negative_inf)
        r = ieee_value(0.0_dp, ieee_positive_inf)
        if (present(lhs)) l = lhs
        if (present(rhs)) r = rhs
        if (present(name)) then
            nm = name
        else
            nc = self%get_nconss()
            write(nm, '("c",i0)') nc + 1
        end if
        cnm = c_string(trim(nm))
        rc = c_scipf_add_linear_cons(self%handle, int(n, c_int), cv, cc, l, r, cnm, idx)
        call return_error(rc, 'scip_add_linear_cons', ierr)
        if (rc == 0_c_int) then
            index1 = int(idx)
        else
            index1 = -1
        end if
    end function

    integer function model_add_quadratic_cons(self, quadvars1, quadvars2, quadcoefs, &
        linvars, lincoefs, lhs, rhs, name, ierr) result(index1)
        class(scip_model_t), intent(inout) :: self
        integer, intent(in) :: quadvars1(:), quadvars2(:)
        real(dp), intent(in) :: quadcoefs(:)
        integer, intent(in), optional :: linvars(:)
        real(dp), intent(in), optional :: lincoefs(:)
        real(dp), intent(in), optional :: lhs, rhs
        character(len=*), intent(in), optional :: name
        integer, intent(out), optional :: ierr
        integer(c_int), allocatable :: lv(:), q1(:), q2(:)
        real(c_double), allocatable :: lc(:), qc(:)
        integer :: nlin, nquad, nc
        real(dp) :: l, r
        character(len=64) :: nm
        character(c_char), allocatable :: cnm(:)
        integer(c_int) :: rc, idx
        call validate_model(self)
        nquad = size(quadvars1)
        if (size(quadvars2) /= nquad .or. size(quadcoefs) /= nquad) &
            error stop 'scip_add_quadratic_cons: quadratic size mismatch'
        if (present(linvars) .neqv. present(lincoefs)) &
            error stop 'scip_add_quadratic_cons: linvars and lincoefs must appear together'
        nlin = 0
        if (present(linvars)) then
            nlin = size(linvars)
            if (size(lincoefs) /= nlin) error stop 'scip_add_quadratic_cons: linear size mismatch'
        end if
        allocate(lv(max(1, nlin)), lc(max(1, nlin)))
        allocate(q1(max(1, nquad)), q2(max(1, nquad)), qc(max(1, nquad)))
        if (nlin > 0) then
            lv(1:nlin) = int(linvars, c_int)
            lc(1:nlin) = real(lincoefs, c_double)
        end if
        if (nquad > 0) then
            q1(1:nquad) = int(quadvars1, c_int)
            q2(1:nquad) = int(quadvars2, c_int)
            qc(1:nquad) = real(quadcoefs, c_double)
        end if
        l = ieee_value(0.0_dp, ieee_negative_inf)
        r = ieee_value(0.0_dp, ieee_positive_inf)
        if (present(lhs)) l = lhs
        if (present(rhs)) r = rhs
        if (present(name)) then
            nm = name
        else
            nc = self%get_nconss()
            write(nm, '("qc",i0)') nc + 1
        end if
        cnm = c_string(trim(nm))
        rc = c_scipf_add_quadratic_cons(self%handle, int(nlin, c_int), lv, lc, &
            int(nquad, c_int), q1, q2, qc, l, r, cnm, idx)
        call return_error(rc, 'scip_add_quadratic_cons', ierr)
        if (rc == 0_c_int) then
            index1 = int(idx)
        else
            index1 = -1
        end if
    end function

    integer function model_add_sos1_cons(self, vars, weights, name, ierr) result(index1)
        class(scip_model_t), intent(inout) :: self
        integer, intent(in) :: vars(:)
        real(dp), intent(in), optional :: weights(:)
        character(len=*), intent(in), optional :: name
        integer, intent(out), optional :: ierr
        index1 = model_add_sos(self, .false., vars, weights, name, ierr)
    end function

    integer function model_add_sos2_cons(self, vars, weights, name, ierr) result(index1)
        class(scip_model_t), intent(inout) :: self
        integer, intent(in) :: vars(:)
        real(dp), intent(in), optional :: weights(:)
        character(len=*), intent(in), optional :: name
        integer, intent(out), optional :: ierr
        index1 = model_add_sos(self, .true., vars, weights, name, ierr)
    end function

    integer function model_add_sos(self, sos2, vars, weights, name, ierr) result(index1)
        class(scip_model_t), intent(inout) :: self
        logical, intent(in) :: sos2
        integer, intent(in) :: vars(:)
        real(dp), intent(in), optional :: weights(:)
        character(len=*), intent(in), optional :: name
        integer, intent(out), optional :: ierr
        integer(c_int), allocatable :: cv(:)
        real(c_double), allocatable :: cw(:)
        integer :: n, j, nc
        integer(c_int) :: rc, idx
        character(len=64) :: nm
        character(c_char), allocatable :: cnm(:)
        call validate_model(self)
        n = size(vars)
        if (present(weights)) then
            if (size(weights) /= n) error stop 'scip_add_sos: size mismatch'
        end if
        allocate(cv(max(1, n)), cw(max(1, n)))
        if (n > 0) cv(1:n) = int(vars, c_int)
        if (present(weights)) then
            if (n > 0) cw(1:n) = real(weights, c_double)
        else
            do j = 1, n
                cw(j) = real(j, c_double)
            end do
        end if
        if (present(name)) then
            nm = name
        else
            nc = self%get_nconss()
            if (sos2) then
                write(nm, '("sos2_",i0)') nc + 1
            else
                write(nm, '("sos1_",i0)') nc + 1
            end if
        end if
        cnm = c_string(trim(nm))
        if (sos2) then
            rc = c_scipf_add_sos2_cons(self%handle, int(n, c_int), cv, cw, cnm, idx)
        else
            rc = c_scipf_add_sos1_cons(self%handle, int(n, c_int), cv, cw, cnm, idx)
        end if
        call return_error(rc, 'scip_add_sos', ierr)
        if (rc == 0_c_int) then
            index1 = int(idx)
        else
            index1 = -1
        end if
    end function

    integer function model_add_indicator_cons(self, binvar, vars, coefs, rhs, name, ierr) &
        result(index1)
        class(scip_model_t), intent(inout) :: self
        integer, intent(in) :: binvar
        integer, intent(in) :: vars(:)
        real(dp), intent(in) :: coefs(:)
        real(dp), intent(in) :: rhs
        character(len=*), intent(in), optional :: name
        integer, intent(out), optional :: ierr
        integer(c_int), allocatable :: cv(:)
        real(c_double), allocatable :: cc(:)
        integer :: n, nc
        integer(c_int) :: rc, idx
        character(len=64) :: nm
        character(c_char), allocatable :: cnm(:)
        call validate_model(self)
        n = size(vars)
        if (size(coefs) /= n) error stop 'scip_add_indicator_cons: size mismatch'
        allocate(cv(max(1, n)), cc(max(1, n)))
        if (n > 0) then
            cv(1:n) = int(vars, c_int)
            cc(1:n) = real(coefs, c_double)
        end if
        if (present(name)) then
            nm = name
        else
            nc = self%get_nconss()
            write(nm, '("ind_",i0)') nc + 1
        end if
        cnm = c_string(trim(nm))
        rc = c_scipf_add_indicator_cons(self%handle, int(binvar, c_int), int(n, c_int), &
            cv, cc, rhs, cnm, idx)
        call return_error(rc, 'scip_add_indicator_cons', ierr)
        if (rc == 0_c_int) then
            index1 = int(idx)
        else
            index1 = -1
        end if
    end function

    subroutine model_set_param_integer(self, name, value, ierr)
        class(scip_model_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer, intent(in) :: value
        integer, intent(out), optional :: ierr
        call model_set_param_int64(self, name, int(value, c_long_long), ierr)
    end subroutine

    subroutine model_set_param_int64(self, name, value, ierr)
        class(scip_model_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer(c_long_long), intent(in) :: value
        integer, intent(out), optional :: ierr
        character(c_char), allocatable :: cname(:)
        integer(c_int) :: rc
        call validate_model(self)
        cname = c_string(name)
        rc = c_scipf_set_param_integer_auto(self%handle, cname, value)
        call return_error(rc, 'scip_set_param(integer)', ierr)
    end subroutine

    subroutine model_set_param_real(self, name, value, ierr)
        class(scip_model_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: value
        integer, intent(out), optional :: ierr
        character(c_char), allocatable :: cname(:)
        integer(c_int) :: rc
        call validate_model(self)
        cname = c_string(name)
        rc = c_scipf_set_param_real_auto(self%handle, cname, value)
        call return_error(rc, 'scip_set_param(real)', ierr)
    end subroutine

    subroutine model_set_param_logical(self, name, value, ierr)
        class(scip_model_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        logical, intent(in) :: value
        integer, intent(out), optional :: ierr
        call model_set_param_int64(self, name, merge(1_c_long_long, 0_c_long_long, value), ierr)
    end subroutine

    subroutine model_set_param_text(self, name, value, ierr)
        class(scip_model_t), intent(inout) :: self
        character(len=*), intent(in) :: name, value
        integer, intent(out), optional :: ierr
        character(c_char), allocatable :: cname(:), cvalue(:)
        integer(c_int) :: rc
        call validate_model(self)
        cname = c_string(name)
        cvalue = c_string(value)
        rc = c_scipf_set_param_text_auto(self%handle, cname, cvalue)
        call return_error(rc, 'scip_set_param(text)', ierr)
    end subroutine

    subroutine model_set_objective_sense(self, sense, ierr)
        class(scip_model_t), intent(inout) :: self
        character(len=*), intent(in) :: sense
        integer, intent(out), optional :: ierr
        integer(c_int) :: rc, imax
        call validate_model(self)
        imax = 0_c_int
        if (sense == 'maximize' .or. sense == 'max') imax = 1_c_int
        rc = c_scipf_set_objective_sense(self%handle, imax)
        call return_error(rc, 'scip_set_objective_sense', ierr)
    end subroutine

    subroutine model_optimize(self, ierr)
        class(scip_model_t), intent(inout) :: self
        integer, intent(out), optional :: ierr
        integer(c_int) :: rc
        call validate_model(self)
        rc = c_scipf_optimize(self%handle)
        call return_error(rc, 'scip_optimize', ierr)
    end subroutine

    function model_get_status(self) result(status)
        class(scip_model_t), intent(in) :: self
        character(len=32) :: status
        character(c_char) :: buf(64)
        integer(c_int) :: rc
        call validate_model(self)
        buf = c_null_char
        rc = c_scipf_get_status(self%handle, buf, int(size(buf), c_int))
        if (rc /= 0_c_int) error stop 'scip_get_status: backend error'
        status = f_string(buf)
    end function

    function model_get_solution(self) result(sol)
        class(scip_model_t), intent(in) :: self
        type(scip_solution) :: sol
        integer :: n
        integer(c_int) :: rc
        real(c_double) :: obj
        real(c_double), allocatable :: x(:)
        call validate_model(self)
        n = self%get_nvars()
        allocate(x(max(1, n)))
        rc = c_scipf_get_best_solution(self%handle, obj, x, int(n, c_int))
        if (rc == 0_c_int) then
            sol%available = .true.
            sol%objval = real(obj, dp)
            allocate(sol%x(n))
            if (n > 0) sol%x = real(x(1:n), dp)
        else
            sol%available = .false.
            sol%objval = ieee_value(0.0_dp, ieee_quiet_nan)
            allocate(sol%x(0))
        end if
    end function

    real(dp) function model_get_objval(self) result(value)
        class(scip_model_t), intent(in) :: self
        type(scip_solution) :: sol
        sol = self%get_solution()
        value = sol%objval
    end function

    integer function model_get_nsols(self) result(n)
        class(scip_model_t), intent(in) :: self
        call validate_model(self)
        n = int(c_scipf_get_nsols(self%handle))
    end function

    function model_get_sol(self, k) result(sol)
        class(scip_model_t), intent(in) :: self
        integer, intent(in) :: k
        type(scip_solution) :: sol
        integer :: n
        integer(c_int) :: rc
        real(c_double) :: obj
        real(c_double), allocatable :: x(:)
        call validate_model(self)
        n = self%get_nvars()
        allocate(x(max(1, n)))
        rc = c_scipf_get_solution_k(self%handle, int(k, c_int), obj, x, int(n, c_int))
        if (rc /= 0_c_int) error stop 'scip_get_sol: solution index out of range'
        sol%available = .true.
        sol%objval = real(obj, dp)
        allocate(sol%x(n))
        if (n > 0) sol%x = real(x(1:n), dp)
    end function

    function model_get_info(self) result(info)
        class(scip_model_t), intent(in) :: self
        type(scip_info) :: info
        real(c_double) :: t, gap
        integer(c_long_long) :: nodes, iters
        integer(c_int) :: ns, rc
        call validate_model(self)
        rc = c_scipf_get_info(self%handle, t, nodes, iters, gap, ns)
        if (rc /= 0_c_int) error stop 'scip_get_info: backend error'
        info%solve_time = real(t, dp)
        info%nodes = nodes
        info%iterations = iters
        info%gap = real(gap, dp)
        info%sol_count = int(ns)
    end function

    integer function model_get_nvars(self) result(n)
        class(scip_model_t), intent(in) :: self
        call validate_model(self)
        n = int(c_scipf_get_nvars(self%handle))
    end function

    integer function model_get_nconss(self) result(n)
        class(scip_model_t), intent(in) :: self
        call validate_model(self)
        n = int(c_scipf_get_nconss(self%handle))
    end function

    subroutine control_append(self, p)
        class(scip_control), intent(inout) :: self
        type(scip_parameter), intent(in) :: p
        type(scip_parameter), allocatable :: tmp(:)
        integer :: n
        n = 0
        if (allocated(self%extra)) n = size(self%extra)
        allocate(tmp(n + 1))
        if (n > 0) tmp(1:n) = self%extra
        tmp(n + 1) = p
        call move_alloc(tmp, self%extra)
    end subroutine

    subroutine control_add_param_logical(self, name, value)
        class(scip_control), intent(inout) :: self
        character(len=*), intent(in) :: name
        logical, intent(in) :: value
        type(scip_parameter) :: p
        p%name = trim(name); p%kind = PARAM_LOGICAL; p%lvalue = value
        call control_append(self, p)
    end subroutine

    subroutine control_add_param_integer(self, name, value)
        class(scip_control), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer, intent(in) :: value
        call control_add_param_int64(self, name, int(value, c_long_long))
    end subroutine

    subroutine control_add_param_int64(self, name, value)
        class(scip_control), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer(c_long_long), intent(in) :: value
        type(scip_parameter) :: p
        p%name = trim(name); p%kind = PARAM_INTEGER; p%ivalue = value
        call control_append(self, p)
    end subroutine

    subroutine control_add_param_real(self, name, value)
        class(scip_control), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: value
        type(scip_parameter) :: p
        p%name = trim(name); p%kind = PARAM_REAL; p%rvalue = value
        call control_append(self, p)
    end subroutine

    subroutine control_add_param_text(self, name, value)
        class(scip_control), intent(inout) :: self
        character(len=*), intent(in) :: name, value
        type(scip_parameter) :: p
        p%name = trim(name); p%kind = PARAM_TEXT; p%tvalue = value
        call control_append(self, p)
    end subroutine

    subroutine apply_control(model, ctrl, ierr)
        type(scip_model_t), intent(inout) :: model
        type(scip_control), intent(in) :: ctrl
        integer, intent(out) :: ierr
        integer :: er, j, hs
        ierr = 0
        if (ctrl%verbose) then
            call model%set_param('display/verblevel', ctrl%verbosity_level, er)
        else
            call model%set_param('display/verblevel', 0, er)
        end if
        if (er /= 0) then; ierr = er; return; end if
        if (ctrl%display_freq /= 100) then
            call model%set_param('display/freq', ctrl%display_freq, er)
            if (er /= 0) then; ierr = er; return; end if
        end if
        if (ctrl%time_limit < huge(1.0_dp) / 2.0_dp) then
            call model%set_param('limits/time', ctrl%time_limit, er)
            if (er /= 0) then; ierr = er; return; end if
        end if
        if (ctrl%node_limit > 0) then
            call model%set_param('limits/nodes', ctrl%node_limit, er)
            if (er /= 0) then; ierr = er; return; end if
        end if
        if (ctrl%stall_node_limit > 0) then
            call model%set_param('limits/stallnodes', ctrl%stall_node_limit, er)
            if (er /= 0) then; ierr = er; return; end if
        end if
        if (ctrl%sol_limit > 0) call model%set_param('limits/solutions', ctrl%sol_limit, er)
        if (er /= 0) then; ierr = er; return; end if
        if (ctrl%best_sol_limit > 0) call model%set_param('limits/bestsol', ctrl%best_sol_limit, er)
        if (er /= 0) then; ierr = er; return; end if
        if (ctrl%mem_limit < huge(1.0_dp) / 2.0_dp) call model%set_param('limits/memory', ctrl%mem_limit, er)
        if (er /= 0) then; ierr = er; return; end if
        if (ctrl%restart_limit > 0) call model%set_param('limits/restarts', ctrl%restart_limit, er)
        if (er /= 0) then; ierr = er; return; end if
        if (ctrl%gap_limit > 0.0_dp) call model%set_param('limits/gap', ctrl%gap_limit, er)
        if (er /= 0) then; ierr = er; return; end if
        if (ctrl%abs_gap_limit > 0.0_dp) call model%set_param('limits/absgap', ctrl%abs_gap_limit, er)
        if (er /= 0) then; ierr = er; return; end if
        if (abs(ctrl%feastol - 1.0e-6_dp) > tiny(1.0_dp)) &
            call model%set_param('numerics/feastol', ctrl%feastol, er)
        if (er /= 0) then; ierr = er; return; end if
        if (abs(ctrl%dualfeastol - 1.0e-7_dp) > tiny(1.0_dp)) &
            call model%set_param('numerics/dualfeastol', ctrl%dualfeastol, er)
        if (er /= 0) then; ierr = er; return; end if
        if (abs(ctrl%epsilon - 1.0e-9_dp) > tiny(1.0_dp)) &
            call model%set_param('numerics/epsilon', ctrl%epsilon, er)
        if (er /= 0) then; ierr = er; return; end if
        if (.not. ctrl%presolving) then
            call model%set_param('presolving/maxrounds', 0, er)
        else if (ctrl%presolve_rounds /= -1) then
            call model%set_param('presolving/maxrounds', ctrl%presolve_rounds, er)
        end if
        if (er /= 0) then; ierr = er; return; end if
        if (ctrl%lp_threads /= 1) call model%set_param('lp/threads', ctrl%lp_threads, er)
        if (er /= 0) then; ierr = er; return; end if
        if (ctrl%lp_iteration_limit > 0) &
            call model%set_param('lp/iterlim', ctrl%lp_iteration_limit, er)
        if (er /= 0) then; ierr = er; return; end if
        if (.not. ctrl%lp_scaling) call model%set_param('lp/scaling', 0, er)
        if (er /= 0) then; ierr = er; return; end if
        if (ctrl%branching_score /= 'p') &
            call model%set_param('branching/scorefunc', ctrl%branching_score, er)
        if (er /= 0) then; ierr = er; return; end if
        if (ctrl%threads /= 1) call model%set_param('parallel/maxnthreads', ctrl%threads, er)
        if (er /= 0) then; ierr = er; return; end if
        hs = 0
        select case (trim(ctrl%heuristics_emphasis))
        case ('aggressive'); hs = 1
        case ('fast'); hs = 2
        case ('off'); hs = 3
        end select
        if (hs /= 0) then
            er = int(c_scipf_set_heuristics(model%handle, int(hs, c_int)))
            if (er /= 0) then; ierr = er; return; end if
        end if
        if (allocated(ctrl%extra)) then
            do j = 1, size(ctrl%extra)
                select case (ctrl%extra(j)%kind)
                case (PARAM_LOGICAL)
                    call model%set_param(ctrl%extra(j)%name, ctrl%extra(j)%lvalue, er)
                case (PARAM_INTEGER)
                    call model%set_param(ctrl%extra(j)%name, ctrl%extra(j)%ivalue, er)
                case (PARAM_REAL)
                    call model%set_param(ctrl%extra(j)%name, ctrl%extra(j)%rvalue, er)
                case (PARAM_TEXT)
                    call model%set_param(ctrl%extra(j)%name, ctrl%extra(j)%tvalue, er)
                end select
                if (er /= 0) then; ierr = er; return; end if
            end do
        end if
    end subroutine

    function make_csc_matrix(a, tol) result(csc)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(in), optional :: tol
        type(scip_csc_matrix) :: csc
        real(dp) :: eps
        integer :: i, j, k, nnz
        eps = 0.0_dp
        if (present(tol)) eps = max(tol, 0.0_dp)
        csc%nrow = size(a, 1)
        csc%ncol = size(a, 2)
        nnz = count(abs(a) > eps)
        allocate(csc%row(nnz), csc%val(nnz), csc%colptr(csc%ncol + 1))
        k = 0
        csc%colptr(1) = 1
        do j = 1, csc%ncol
            do i = 1, csc%nrow
                if (abs(a(i,j)) > eps) then
                    k = k + 1
                    csc%row(k) = i
                    csc%val(k) = a(i,j)
                end if
            end do
            csc%colptr(j + 1) = k + 1
        end do
    end function

    function scip_solve_dense(obj, a, b, sense, vtype, lb, ub, control) result(res)
        real(dp), intent(in) :: obj(:), a(:,:), b(:)
        character(len=*), intent(in) :: sense(:)
        character(len=*), intent(in), optional :: vtype(:)
        real(dp), intent(in), optional :: lb(:), ub(:)
        type(scip_control), intent(in), optional :: control
        type(scip_result) :: res
        type(scip_csc_matrix) :: csc
        csc = make_csc_matrix(a)
        res = scip_solve_csc(obj, csc, b, sense, vtype, lb, ub, control)
    end function

    function scip_solve_csc(obj, a, b, sense, vtype, lb, ub, control) result(res)
        real(dp), intent(in) :: obj(:), b(:)
        type(scip_csc_matrix), intent(in) :: a
        character(len=*), intent(in) :: sense(:)
        character(len=*), intent(in), optional :: vtype(:)
        real(dp), intent(in), optional :: lb(:), ub(:)
        type(scip_control), intent(in), optional :: control
        type(scip_result) :: res
        type(scip_model_t) :: model
        type(scip_control) :: ctrl
        type(scip_solution) :: sol
        integer, allocatable :: rowcnt(:), rowptr(:), pos(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: i, j, k, p, n, m, er, idx
        real(dp) :: lhs, rhs, l, u
        character(len=1) :: vt
        character(len=64) :: nm
        n = size(obj)
        m = a%nrow
        if (a%ncol /= n) error stop 'scip_solve: A column count does not match obj'
        if (size(b) /= m .or. size(sense) /= m) error stop 'scip_solve: constraint size mismatch'
        if (present(vtype)) then
            if (size(vtype) /= n) error stop 'scip_solve: vtype size mismatch'
        end if
        if (present(lb)) then
            if (size(lb) /= n) error stop 'scip_solve: lb size mismatch'
        end if
        if (present(ub)) then
            if (size(ub) /= n) error stop 'scip_solve: ub size mismatch'
        end if
        if (size(a%colptr) /= n + 1) error stop 'scip_solve: invalid CSC colptr'
        if (size(a%row) /= size(a%val)) error stop 'scip_solve: invalid CSC arrays'
        model = scip_model('scip_fortran', er)
        if (er /= 0) error stop 'scip_solve: model creation failed'
        ctrl = scip_control()
        if (present(control)) ctrl = control
        call apply_control(model, ctrl, er)
        if (er /= 0) error stop 'scip_solve: invalid control parameter'
        do j = 1, n
            l = 0.0_dp
            u = ieee_value(0.0_dp, ieee_positive_inf)
            vt = 'C'
            if (present(lb)) l = lb(j)
            if (present(ub)) u = ub(j)
            if (present(vtype)) vt = vtype(j)(1:1)
            write(nm, '("x",i0)') j - 1
            idx = model%add_var(obj(j), l, u, vt, trim(nm), er)
            if (er /= 0) error stop 'scip_solve: variable creation failed'
        end do
        allocate(rowcnt(m), rowptr(m + 1), pos(m))
        rowcnt = 0
        do k = 1, size(a%row)
            if (a%row(k) < 1 .or. a%row(k) > m) error stop 'scip_solve: CSC row out of range'
            rowcnt(a%row(k)) = rowcnt(a%row(k)) + 1
        end do
        rowptr(1) = 1
        do i = 1, m
            rowptr(i + 1) = rowptr(i) + rowcnt(i)
        end do
        allocate(cols(size(a%row)), vals(size(a%val)))
        pos = rowptr(1:m)
        do j = 1, n
            if (a%colptr(j) < 1 .or. a%colptr(j + 1) < a%colptr(j)) &
                error stop 'scip_solve: invalid CSC colptr'
            do k = a%colptr(j), a%colptr(j + 1) - 1
                i = a%row(k)
                p = pos(i)
                cols(p) = j
                vals(p) = a%val(k)
                pos(i) = p + 1
            end do
        end do
        do i = 1, m
            lhs = ieee_value(0.0_dp, ieee_negative_inf)
            rhs = ieee_value(0.0_dp, ieee_positive_inf)
            select case (trim(sense(i)))
            case ('<=', '<')
                rhs = b(i)
            case ('>=', '>')
                lhs = b(i)
            case ('==', '=')
                lhs = b(i); rhs = b(i)
            case default
                error stop 'scip_solve: sense must be <=, >=, or =='
            end select
            write(nm, '("c",i0)') i - 1
            idx = model%add_linear_cons(cols(rowptr(i):rowptr(i + 1) - 1), &
                vals(rowptr(i):rowptr(i + 1) - 1), lhs, rhs, trim(nm), er)
            if (er /= 0) error stop 'scip_solve: constraint creation failed'
        end do
        call model%optimize(er)
        if (er /= 0) error stop 'scip_solve: optimization failed'
        res%status = model%get_status()
        sol = model%get_solution()
        res%sol_count = model%get_nsols()
        res%info = model%get_info()
        res%gap = res%info%gap
        if (sol%available) then
            res%objval = sol%objval
            allocate(res%x(size(sol%x)))
            res%x = sol%x
        else
            res%objval = ieee_value(0.0_dp, ieee_quiet_nan)
            allocate(res%x(0))
        end if
        call model%free()
    end function

end module scip
