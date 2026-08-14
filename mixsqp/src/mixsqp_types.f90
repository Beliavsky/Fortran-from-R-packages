module mixsqp_types
  use mixsqp_kinds, only : dp
  implicit none
  private
  public :: mixsqp_control, mixsqp_result

  type :: mixsqp_control
    logical :: normalize_rows = .true.
    real(dp) :: tol_svd = 1.0e-6_dp
    real(dp) :: convtol_sqp = 1.0e-8_dp
    real(dp) :: convtol_activeset = 1.0e-10_dp
    real(dp) :: zero_threshold_solution = 1.0e-8_dp
    real(dp) :: zero_threshold_searchdir = 1.0e-14_dp
    real(dp) :: suffdecr_linesearch = 1.0e-2_dp
    real(dp) :: stepsize_reduce = 0.75_dp
    real(dp) :: min_stepsize = 1.0e-8_dp
    real(dp) :: identity_contrib_increase = 10.0_dp
    real(dp) :: eps = 1.0e-8_dp
    integer :: maxiter_sqp = 1000
    integer :: maxiter_activeset = 0
    integer :: numiter_em = 10
    logical :: verbose = .false.
  end type mixsqp_control

  type :: mixsqp_result
    real(dp), allocatable :: x(:)
    real(dp) :: value = huge(1.0_dp)
    real(dp), allocatable :: grad(:)
    real(dp), allocatable :: hessian(:,:)
    integer :: status = 1
    character(len=48) :: status_message = "exceeded maximum number of iterations"
    integer :: iterations = 0
    logical :: used_svd = .false.
    integer :: svd_rank = 0
    real(dp), allocatable :: objective(:)
    real(dp), allocatable :: max_rdual(:)
    integer, allocatable :: nnz(:)
    real(dp), allocatable :: stepsize(:)
    real(dp), allocatable :: max_diff(:)
    integer, allocatable :: nqp(:)
    integer, allocatable :: nls(:)
  end type mixsqp_result
end module mixsqp_types
