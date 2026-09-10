! SPDX-License-Identifier: GPL-2.0-or-later
module dlm_mle_mod
    use dlm_types, only : dp, dlm_model, dlm_success, dlm_invalid_shape, dlm_invalid_argument
    use dlm_kalman, only : dlm_ll
    use optimx_mod, only : optimx_problem, optimx_control, optimx_result
    use optimx_mod, only : initialize_problem, optimr, OPTIMX_INVALID_INPUT
    implicit none
    private

    type, abstract, public :: dlm_model_builder
    contains
        procedure(dlm_model_build_interface), deferred, pass :: build
    end type dlm_model_builder

    type, public :: dlm_mle_result
        real(dp), allocatable :: par(:)
        real(dp) :: value = huge(1.0_dp)
        integer :: function_count = 0
        integer :: gradient_count = 0
        integer :: iterations = 0
        integer :: convergence = OPTIMX_INVALID_INPUT
        logical :: converged = .false.
        character(len=32) :: method = ''
        character(len=160) :: message = 'not run'
    end type dlm_mle_result

    abstract interface
        subroutine dlm_model_build_interface(self, parameters, model, info)
            import :: dlm_model_builder, dlm_model, dp
            class(dlm_model_builder), intent(in) :: self !! Builder object holding any fixed model-construction data.
            real(dp), intent(in) :: parameters(:) !! Free parameter vector used to construct one candidate DLM.
            type(dlm_model), intent(out) :: model !! Candidate dynamic linear model built from `parameters`.
            integer, intent(out) :: info !! Zero when the candidate model is valid or a nonzero builder status otherwise.
        end subroutine dlm_model_build_interface
    end interface

    class(dlm_model_builder), pointer, save :: active_builder => null()
    real(dp), pointer, save :: active_y(:, :) => null()

    interface dlm_mle
        module procedure dlm_mle_matrix
        module procedure dlm_mle_vector
    end interface dlm_mle

    public :: dlm_mle

contains

    subroutine dlm_mle_matrix(y, initial, builder, result, info, method, lower, upper, &
                              max_iterations, max_evaluations, tolerance)
        real(dp), intent(in), target :: y(:, :) !! Observation matrix with time along rows; IEEE NaNs mark missing components.
        real(dp), intent(in) :: initial(:) !! Starting parameter vector passed first to the model builder.
        class(dlm_model_builder), intent(in), target :: builder !! Typed DLM builder used for each parameter vector.
        type(dlm_mle_result), intent(out) :: result !! Optimizer parameters, likelihood value, counts, and convergence metadata.
        integer, intent(out) :: info !! Zero when optimization was launched successfully or a `dlm_*` input status otherwise.
        character(len=*), intent(in), optional :: method !! Optimizer name; defaults to `L-BFGS-B` as in upstream `dlmMLE`.
        real(dp), intent(in), optional :: lower(:) !! Optional lower bounds with one entry per free parameter.
        real(dp), intent(in), optional :: upper(:) !! Optional upper bounds with one entry per free parameter.
        integer, intent(in), optional :: max_iterations !! Optional maximum optimizer iterations; must be positive.
        integer, intent(in), optional :: max_evaluations !! Optional maximum objective evaluations; must be positive.
        real(dp), intent(in), optional :: tolerance !! Optional relative convergence tolerance; must be positive.

        type(optimx_problem) :: problem
        type(optimx_control) :: control
        type(optimx_result) :: optimizer_result
        character(len=32) :: method_name
        integer :: n_parameters

        result = dlm_mle_result()
        info = dlm_success
        n_parameters = size(initial)
        if (n_parameters < 1 .or. size(y, 1) < 1 .or. size(y, 2) < 1) then
            info = dlm_invalid_argument
            return
        end if
        if (present(lower)) then
            if (size(lower) /= n_parameters) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(upper)) then
            if (size(upper) /= n_parameters) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(max_iterations)) then
            if (max_iterations < 1) then
                info = dlm_invalid_argument
                return
            end if
        end if
        if (present(max_evaluations)) then
            if (max_evaluations < 1) then
                info = dlm_invalid_argument
                return
            end if
        end if
        if (present(tolerance)) then
            if (tolerance <= 0.0_dp) then
                info = dlm_invalid_argument
                return
            end if
        end if

        call initialize_problem(problem, n_parameters)
        problem%objective => likelihood_callback
        if (present(lower)) problem%lower = lower
        if (present(upper)) problem%upper = upper
        if (any(problem%lower > problem%upper)) then
            info = dlm_invalid_argument
            return
        end if

        control = optimx_control()
        if (present(max_iterations)) control%maxit = max_iterations
        if (present(max_evaluations)) control%maxfeval = max_evaluations
        if (present(tolerance)) control%reltol = tolerance
        method_name = 'L-BFGS-B'
        if (present(method)) method_name = method
        if (associated(active_builder) .or. associated(active_y)) then
            info = dlm_invalid_argument
            return
        end if

        active_builder => builder
        active_y => y
        call optimr(problem, initial, trim(method_name), control, optimizer_result)
        nullify(active_builder)
        nullify(active_y)
        if (allocated(optimizer_result%par)) result%par = optimizer_result%par
        result%value = optimizer_result%value
        result%function_count = optimizer_result%function_count
        result%gradient_count = optimizer_result%gradient_count
        result%iterations = optimizer_result%iterations
        result%convergence = optimizer_result%convergence
        result%converged = optimizer_result%converged
        result%method = optimizer_result%method
        result%message = optimizer_result%message
        if (optimizer_result%convergence == OPTIMX_INVALID_INPUT) info = dlm_invalid_argument


    end subroutine dlm_mle_matrix

    subroutine likelihood_callback(parameters, value, gradient, hessian, need_gradient, need_hessian, status)
        real(dp), intent(in) :: parameters(:) !! Candidate optimizer parameter vector.
        real(dp), intent(out) :: value !! Constant-free Gaussian negative log likelihood returned to the optimizer.
        real(dp), intent(inout) :: gradient(:) !! Gradient workspace; zeroed because derivatives are generated numerically.
        real(dp), intent(inout) :: hessian(:, :) !! Hessian workspace; zeroed because no analytic Hessian is supplied.
        logical, intent(in) :: need_gradient !! Optimizer request flag for an analytic gradient.
        logical, intent(in) :: need_hessian !! Optimizer request flag for an analytic Hessian.
        integer, intent(out) :: status !! Zero after a valid model/likelihood evaluation or a nonzero DLM status otherwise.

        type(dlm_model) :: model

        gradient = 0.0_dp
        hessian = 0.0_dp
        if (need_gradient .or. need_hessian) continue
        if (.not. associated(active_builder) .or. .not. associated(active_y)) then
            status = dlm_invalid_argument
            value = huge(1.0_dp)
            return
        end if
        call active_builder%build(parameters, model, status)
        if (status /= dlm_success) then
            value = huge(1.0_dp)
            return
        end if
        call dlm_ll(active_y, model, value, status)
        if (status /= dlm_success) value = huge(1.0_dp)
    end subroutine likelihood_callback

    subroutine dlm_mle_vector(y, initial, builder, result, info, method, lower, upper, &
                              max_iterations, max_evaluations, tolerance)
        real(dp), intent(in) :: y(:) !! Univariate observation sequence; IEEE NaNs mark missing observations.
        real(dp), intent(in) :: initial(:) !! Starting parameter vector passed first to the model builder.
        class(dlm_model_builder), intent(in), target :: builder !! Typed DLM builder used for each parameter vector.
        type(dlm_mle_result), intent(out) :: result !! Optimizer parameters, likelihood value, counts, and convergence metadata.
        integer, intent(out) :: info !! Zero when optimization was launched successfully or a `dlm_*` input status otherwise.
        character(len=*), intent(in), optional :: method !! Optimizer name; defaults to `L-BFGS-B` as in upstream `dlmMLE`.
        real(dp), intent(in), optional :: lower(:) !! Optional lower bounds with one entry per free parameter.
        real(dp), intent(in), optional :: upper(:) !! Optional upper bounds with one entry per free parameter.
        integer, intent(in), optional :: max_iterations !! Optional maximum optimizer iterations; must be positive.
        integer, intent(in), optional :: max_evaluations !! Optional maximum objective evaluations; must be positive.
        real(dp), intent(in), optional :: tolerance !! Optional relative convergence tolerance; must be positive.

        real(dp), allocatable, target :: y_matrix(:, :)

        allocate(y_matrix(size(y), 1))
        y_matrix(:, 1) = y
        call dlm_mle_matrix(y_matrix, initial, builder, result, info, method, lower, upper, &
                            max_iterations, max_evaluations, tolerance)
    end subroutine dlm_mle_vector

end module dlm_mle_mod
