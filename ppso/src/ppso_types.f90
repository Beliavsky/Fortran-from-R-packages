module ppso_types
    use ppso_kinds, only : dp
    use ppso_rng, only : rng_state
    implicit none
    private

    integer, parameter, public :: PPSO_STOP_NONE = 0
    integer, parameter, public :: PPSO_STOP_MAX_CALLS = 1
    integer, parameter, public :: PPSO_STOP_MAX_ITER = 2
    integer, parameter, public :: PPSO_STOP_CONVERGED = 3
    integer, parameter, public :: PPSO_STOP_INVALID_OBJECTIVE = 4

    type, public :: pso_control
        integer :: number_of_particles = 40
        integer :: max_number_of_iterations = 5
        integer :: max_number_function_calls = 500
        integer :: max_wait_iterations = 50
        real(dp) :: w = 1.0_dp
        real(dp) :: c1 = 2.0_dp
        real(dp) :: c2 = 2.0_dp
        real(dp) :: abstol = -huge(1.0_dp)
        real(dp) :: reltol = -huge(1.0_dp)
        logical :: wait_complete_iteration = .false.
        logical :: lhs_init = .false.
        integer(kind=8) :: seed = 123456789_8
    end type pso_control

    type, public :: dds_control
        integer :: number_of_particles = 1
        integer :: max_number_function_calls = 500
        integer :: max_wait_iterations = 50
        real(dp) :: r = 0.2_dp
        real(dp) :: abstol = -huge(1.0_dp)
        real(dp) :: reltol = -huge(1.0_dp)
        integer :: part_xchange = 2
        logical :: lhs_init = .false.
        logical :: legacy_serial_prerun_omission = .true.
        integer(kind=8) :: seed = 123456789_8
    end type dds_control

    type, public :: ppso_result
        real(dp) :: value = huge(1.0_dp)
        real(dp), allocatable :: par(:)
        integer :: function_calls = 0
        integer :: actual_function_calls = 0
        integer :: iterations = 0
        integer :: stop_code = PPSO_STOP_NONE
        character(len=80) :: break_flag = ""
    end type ppso_result

    type, public :: pso_state
        integer :: npar = 0
        integer :: npart = 0
        real(dp), allocatable :: lower(:), upper(:), vmax(:)
        real(dp), allocatable :: x(:,:), v(:,:), fitness_x(:)
        real(dp), allocatable :: x_lbest(:,:), fitness_lbest(:)
        real(dp), allocatable :: x_gbest(:)
        real(dp), allocatable :: pending_initial(:,:)
        integer, allocatable :: function_calls(:)
        logical, allocatable :: fixed_until_evaluated(:)
        real(dp) :: fitness_gbest = huge(1.0_dp)
        real(dp) :: fitness_itbest = huge(1.0_dp)
        integer :: it_last_improvement = 0
        integer :: stop_code = PPSO_STOP_NONE
        type(rng_state) :: rng
    end type pso_state

    type, public :: dds_state
        integer :: npar = 0
        integer :: npart = 0
        real(dp), allocatable :: lower(:), upper(:)
        real(dp), allocatable :: x(:,:), v(:,:), fitness_x(:)
        real(dp), allocatable :: x_lbest(:,:), fitness_lbest(:)
        real(dp), allocatable :: x_gbest(:)
        real(dp), allocatable :: pending_initial(:,:)
        integer, allocatable :: function_calls(:), futile_iter_count(:)
        real(dp) :: fitness_gbest = huge(1.0_dp)
        real(dp) :: fitness_itbest = huge(1.0_dp)
        integer :: it_last_improvement = 0
        integer :: init_function_calls = 0
        integer :: actual_init_calls = 0
        integer :: main_max_calls = 0
        integer :: stop_code = PPSO_STOP_NONE
        type(rng_state) :: rng
    end type dds_state

    abstract interface
        function objective_function(x) result(value)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function objective_function

        subroutine ppso_monitor(iteration, calls, best_value, best_par)
            import dp
            integer, intent(in) :: iteration, calls
            real(dp), intent(in) :: best_value
            real(dp), intent(in) :: best_par(:)
        end subroutine ppso_monitor
    end interface

    public :: objective_function, ppso_monitor, stop_message

contains

    function stop_message(code) result(msg)
        integer, intent(in) :: code
        character(len=80) :: msg

        select case (code)
        case (PPSO_STOP_MAX_CALLS)
            msg = "max number of function calls reached"
        case (PPSO_STOP_MAX_ITER)
            msg = "max iterations reached"
        case (PPSO_STOP_CONVERGED)
            msg = "converged"
        case (PPSO_STOP_INVALID_OBJECTIVE)
            msg = "objective function returned NaN"
        case default
            msg = ""
        end select
    end function stop_message

end module ppso_types
