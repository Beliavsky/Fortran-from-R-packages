module deoptim_types
    use deoptim_kinds, only : dp, i8
    implicit none
    private

    integer, parameter, public :: de_success = 0
    integer, parameter, public :: de_invalid_input = 1
    integer, parameter, public :: de_unsupported = 2
    integer, parameter, public :: de_objective_nan = 3
    integer, parameter, public :: de_map_error = 4

    type, public :: de_control
        real(dp) :: vtr = -huge(1.0_dp)
        integer :: strategy = 2
        integer :: np = 0
        integer :: itermax = 200
        real(dp) :: cr = 0.5_dp
        real(dp) :: f = 0.8_dp
        logical :: bs = .false.
        integer :: trace = 0
        integer :: storepopfrom = 0
        integer :: storepopfreq = 1
        real(dp) :: p = 0.2_dp
        real(dp) :: c = 0.0_dp
        real(dp) :: reltol = sqrt(epsilon(1.0_dp))
        integer :: steptol = 0
        integer(i8) :: seed = 0_i8
    end type de_control

    type, public :: de_result
        real(dp), allocatable :: bestmem(:)
        real(dp) :: bestval = huge(1.0_dp)
        integer(i8) :: nfeval = 0_i8
        integer :: iter = 0
        real(dp), allocatable :: bestmemit(:,:)
        real(dp), allocatable :: bestvalit(:)
        real(dp), allocatable :: pop(:,:)
        real(dp), allocatable :: storepop(:,:,:)
        integer :: nstore = 0
        integer :: status = de_success
        character(len=:), allocatable :: message
        integer :: map_duplicates_remaining = 0
    end type de_result

    abstract interface
        function de_objective(x) result(value)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function de_objective

        subroutine de_map(x)
            import dp
            real(dp), intent(inout) :: x(:)
        end subroutine de_map
    end interface

    public :: de_objective, de_map
end module deoptim_types
