! SPDX-License-Identifier: GPL-3.0-only
module qpoases
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
    use qpoases_kinds, only : dp
    use qpoases_types
    use qpoases_solver, only : solve_general_qp
    implicit none
    private

    public :: dp
    public :: qpoases_options, qpoases_result, qpoases_model
    public :: hst_zero, hst_identity, hst_posdef, hst_posdef_nullspace
    public :: hst_semidef, hst_indef, hst_unknown
    public :: successful_return, ret_invalid_arguments, ret_qp_unbounded
    public :: ret_qp_infeasible, ret_qp_not_solved, ret_qp_solved
    public :: ret_hotstart_failed, ret_max_nwsr_reached
    public :: solve_qproblem, solve_qproblemb
    public :: init_qproblem, init_qproblemb, init_sqproblem
    public :: hotstart_qproblem, hotstart_qproblemb, hotstart_sqproblem
    public :: set_options, default_options
    public :: get_objval, get_primal_solution, get_dual_solution
    public :: get_number_of_variables, get_number_of_free_variables
    public :: get_number_of_fixed_variables, get_number_of_constraints
    public :: get_number_of_equality_constraints
    public :: get_number_of_active_constraints, get_number_of_inactive_constraints
    public :: is_initialised, is_solved, is_infeasible, is_unbounded
    public :: status_symbol

    type, public :: qpoases_model
        integer :: n_variables = 0
        integer :: n_constraints = 0
        integer :: hessian_type = hst_unknown
        real(dp), allocatable :: h(:,:), g(:), a(:,:)
        real(dp), allocatable :: lb(:), ub(:), lba(:), uba(:)
        type(qpoases_options) :: options
        type(qpoases_result) :: result
        integer, allocatable :: active_ids(:)
    end type qpoases_model

contains

    pure function default_options() result(options)
        type(qpoases_options) :: options
        options = qpoases_options()
    end function default_options

    subroutine set_options(model, options)
        type(qpoases_model), intent(inout) :: model
        type(qpoases_options), intent(in) :: options
        model%options = options
    end subroutine set_options

    subroutine solve_qproblem(h, g, a, lb, ub, lba, uba, result, options, &
                              hessian_type, max_nwsr)
        real(dp), intent(in) :: h(:,:), g(:), a(:,:), lb(:), ub(:), lba(:), uba(:)
        type(qpoases_result), intent(out) :: result
        type(qpoases_options), intent(in), optional :: options
        integer, intent(in), optional :: hessian_type, max_nwsr
        type(qpoases_options) :: opt
        integer :: ht, nwsr

        opt = default_options()
        if (present(options)) opt = options
        ht = hst_unknown
        if (present(hessian_type)) ht = hessian_type
        nwsr = 2000
        if (present(max_nwsr)) nwsr = max_nwsr
        call solve_general_qp(h,g,a,lb,ub,lba,uba,ht,opt,nwsr,result)
    end subroutine solve_qproblem

    subroutine solve_qproblemb(h, g, lb, ub, result, options, hessian_type, max_nwsr)
        real(dp), intent(in) :: h(:,:), g(:), lb(:), ub(:)
        type(qpoases_result), intent(out) :: result
        type(qpoases_options), intent(in), optional :: options
        integer, intent(in), optional :: hessian_type, max_nwsr
        real(dp), allocatable :: a(:,:), lba(:), uba(:)
        type(qpoases_options) :: opt
        integer :: ht, nwsr, n

        n = size(g)
        allocate(a(0,n),lba(0),uba(0))
        opt = default_options()
        if (present(options)) opt = options
        ht = hst_unknown
        if (present(hessian_type)) ht = hessian_type
        nwsr = 2000
        if (present(max_nwsr)) nwsr = max_nwsr
        call solve_general_qp(h,g,a,lb,ub,lba,uba,ht,opt,nwsr,result)
    end subroutine solve_qproblemb

    subroutine store_problem(model, h, g, a, lb, ub, lba, uba, hessian_type)
        type(qpoases_model), intent(inout) :: model
        real(dp), intent(in) :: h(:,:), g(:), a(:,:), lb(:), ub(:), lba(:), uba(:)
        integer, intent(in) :: hessian_type

        model%n_variables = size(g)
        model%n_constraints = size(lba)
        model%hessian_type = hessian_type
        model%h = h
        model%g = g
        model%a = a
        model%lb = lb
        model%ub = ub
        model%lba = lba
        model%uba = uba
    end subroutine store_problem

    subroutine init_qproblem(model, h, g, a, lb, ub, lba, uba, max_nwsr, &
                             status, hessian_type, options)
        type(qpoases_model), intent(inout) :: model
        real(dp), intent(in) :: h(:,:), g(:), a(:,:), lb(:), ub(:), lba(:), uba(:)
        integer, intent(in) :: max_nwsr
        integer, intent(out) :: status
        integer, intent(in), optional :: hessian_type
        type(qpoases_options), intent(in), optional :: options
        integer :: ht

        ht = hst_unknown
        if (present(hessian_type)) ht = hessian_type
        if (present(options)) model%options = options
        call store_problem(model,h,g,a,lb,ub,lba,uba,ht)
        call solve_general_qp(model%h,model%g,model%a,model%lb,model%ub, &
            model%lba,model%uba,model%hessian_type,model%options,max_nwsr, &
            model%result,active_ids=model%active_ids)
        status = model%result%status
    end subroutine init_qproblem

    subroutine init_qproblemb(model, h, g, lb, ub, max_nwsr, status, &
                              hessian_type, options)
        type(qpoases_model), intent(inout) :: model
        real(dp), intent(in) :: h(:,:), g(:), lb(:), ub(:)
        integer, intent(in) :: max_nwsr
        integer, intent(out) :: status
        integer, intent(in), optional :: hessian_type
        type(qpoases_options), intent(in), optional :: options
        real(dp), allocatable :: a(:,:), lba(:), uba(:)
        integer :: n

        n = size(g)
        allocate(a(0,n),lba(0),uba(0))
        call init_qproblem(model,h,g,a,lb,ub,lba,uba,max_nwsr,status, &
                           hessian_type,options)
    end subroutine init_qproblemb

    subroutine init_sqproblem(model, h, g, a, lb, ub, lba, uba, max_nwsr, &
                              status, hessian_type, options)
        type(qpoases_model), intent(inout) :: model
        real(dp), intent(in) :: h(:,:), g(:), a(:,:), lb(:), ub(:), lba(:), uba(:)
        integer, intent(in) :: max_nwsr
        integer, intent(out) :: status
        integer, intent(in), optional :: hessian_type
        type(qpoases_options), intent(in), optional :: options
        call init_qproblem(model,h,g,a,lb,ub,lba,uba,max_nwsr,status, &
                           hessian_type,options)
    end subroutine init_sqproblem

    subroutine hotstart_qproblem(model, g, lb, ub, lba, uba, max_nwsr, status)
        type(qpoases_model), intent(inout) :: model
        real(dp), intent(in) :: g(:), lb(:), ub(:), lba(:), uba(:)
        integer, intent(in) :: max_nwsr
        integer, intent(out) :: status
        real(dp), allocatable :: xstart(:)
        integer, allocatable :: warm(:)

        if (.not. model%result%initialised) then
            status = ret_hotstart_failed
            return
        end if
        if (size(g) /= model%n_variables .or. size(lb) /= model%n_variables .or. &
            size(ub) /= model%n_variables .or. size(lba) /= model%n_constraints .or. &
            size(uba) /= model%n_constraints) then
            status = ret_invalid_arguments
            return
        end if
        xstart = model%result%x
        if (allocated(model%active_ids)) warm = model%active_ids
        model%g = g
        model%lb = lb
        model%ub = ub
        model%lba = lba
        model%uba = uba
        if (allocated(warm)) then
            call solve_general_qp(model%h,model%g,model%a,model%lb,model%ub, &
                model%lba,model%uba,model%hessian_type,model%options,max_nwsr, &
                model%result,xstart,warm,model%active_ids)
        else
            call solve_general_qp(model%h,model%g,model%a,model%lb,model%ub, &
                model%lba,model%uba,model%hessian_type,model%options,max_nwsr, &
                model%result,xstart,active_ids=model%active_ids)
        end if
        status = model%result%status
    end subroutine hotstart_qproblem

    subroutine hotstart_qproblemb(model, g, lb, ub, max_nwsr, status)
        type(qpoases_model), intent(inout) :: model
        real(dp), intent(in) :: g(:), lb(:), ub(:)
        integer, intent(in) :: max_nwsr
        integer, intent(out) :: status
        real(dp), allocatable :: lba(:), uba(:)

        allocate(lba(0),uba(0))
        call hotstart_qproblem(model,g,lb,ub,lba,uba,max_nwsr,status)
    end subroutine hotstart_qproblemb

    subroutine hotstart_sqproblem(model, h, g, a, lb, ub, lba, uba, &
                                  max_nwsr, status)
        type(qpoases_model), intent(inout) :: model
        real(dp), intent(in) :: h(:,:), g(:), a(:,:), lb(:), ub(:), lba(:), uba(:)
        integer, intent(in) :: max_nwsr
        integer, intent(out) :: status
        real(dp), allocatable :: xstart(:)
        integer, allocatable :: warm(:)

        if (.not. model%result%initialised) then
            status = ret_hotstart_failed
            return
        end if
        xstart = model%result%x
        if (allocated(model%active_ids)) warm = model%active_ids
        call store_problem(model,h,g,a,lb,ub,lba,uba,model%hessian_type)
        if (allocated(warm)) then
            call solve_general_qp(model%h,model%g,model%a,model%lb,model%ub, &
                model%lba,model%uba,model%hessian_type,model%options,max_nwsr, &
                model%result,xstart,warm,model%active_ids)
        else
            call solve_general_qp(model%h,model%g,model%a,model%lb,model%ub, &
                model%lba,model%uba,model%hessian_type,model%options,max_nwsr, &
                model%result,xstart,active_ids=model%active_ids)
        end if
        status = model%result%status
    end subroutine hotstart_sqproblem

    real(dp) function get_objval(model) result(value)
        type(qpoases_model), intent(in) :: model
        value = model%result%objval
    end function get_objval

    subroutine get_primal_solution(model, x)
        type(qpoases_model), intent(in) :: model
        real(dp), allocatable, intent(out) :: x(:)
        x = model%result%x
    end subroutine get_primal_solution

    subroutine get_dual_solution(model, y)
        type(qpoases_model), intent(in) :: model
        real(dp), allocatable, intent(out) :: y(:)
        y = model%result%y
    end subroutine get_dual_solution

    integer function get_number_of_variables(model) result(n)
        type(qpoases_model), intent(in) :: model
        n = model%n_variables
    end function get_number_of_variables

    integer function get_number_of_free_variables(model) result(n)
        type(qpoases_model), intent(in) :: model
        n = model%result%n_free
    end function get_number_of_free_variables

    integer function get_number_of_fixed_variables(model) result(n)
        type(qpoases_model), intent(in) :: model
        n = model%result%n_fixed
    end function get_number_of_fixed_variables

    integer function get_number_of_constraints(model) result(n)
        type(qpoases_model), intent(in) :: model
        n = model%n_constraints
    end function get_number_of_constraints

    integer function get_number_of_equality_constraints(model) result(n)
        type(qpoases_model), intent(in) :: model
        integer :: i
        real(dp) :: tol, scale
        n = 0
        tol = max(model%options%bound_tolerance,100.0_dp*epsilon(1.0_dp))
        do i = 1, model%n_constraints
            if (.not. ieee_is_finite_local(model%lba(i))) cycle
            if (.not. ieee_is_finite_local(model%uba(i))) cycle
            scale = max(1.0_dp,abs(model%lba(i)),abs(model%uba(i)))
            if (abs(model%uba(i)-model%lba(i)) <= tol*scale) n = n + 1
        end do
    end function get_number_of_equality_constraints

    pure logical function ieee_is_finite_local(x)
        use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
        real(dp), intent(in) :: x
        ieee_is_finite_local = ieee_is_finite(x)
    end function ieee_is_finite_local

    integer function get_number_of_active_constraints(model) result(n)
        type(qpoases_model), intent(in) :: model
        n = model%result%n_active_constraints
    end function get_number_of_active_constraints

    integer function get_number_of_inactive_constraints(model) result(n)
        type(qpoases_model), intent(in) :: model
        n = model%result%n_inactive_constraints
    end function get_number_of_inactive_constraints

    logical function is_initialised(model)
        type(qpoases_model), intent(in) :: model
        is_initialised = model%result%initialised
    end function is_initialised

    logical function is_solved(model)
        type(qpoases_model), intent(in) :: model
        is_solved = model%result%solved
    end function is_solved

    logical function is_infeasible(model)
        type(qpoases_model), intent(in) :: model
        is_infeasible = model%result%infeasible
    end function is_infeasible

    logical function is_unbounded(model)
        type(qpoases_model), intent(in) :: model
        is_unbounded = model%result%unbounded
    end function is_unbounded

    function status_symbol(status) result(symbol)
        integer, intent(in) :: status
        character(len=:), allocatable :: symbol
        select case (status)
        case (successful_return)
            symbol = "SUCCESSFUL_RETURN"
        case (ret_invalid_arguments)
            symbol = "RET_INVALID_ARGUMENTS"
        case (ret_qp_unbounded)
            symbol = "RET_QP_UNBOUNDED"
        case (ret_qp_infeasible)
            symbol = "RET_QP_INFEASIBLE"
        case (ret_qp_not_solved)
            symbol = "RET_QP_NOT_SOLVED"
        case (ret_qp_solved)
            symbol = "RET_QP_SOLVED"
        case (ret_hotstart_failed)
            symbol = "RET_HOTSTART_FAILED"
        case (ret_max_nwsr_reached)
            symbol = "RET_MAX_NWSR_REACHED"
        case default
            symbol = "RET_STATUS_" // int_to_string(status)
        end select
    end function status_symbol

    function int_to_string(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=32) :: buffer
        write(buffer,'(i0)') i
        s = trim(buffer)
    end function int_to_string
end module qpoases
