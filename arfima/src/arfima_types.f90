module arfima_types
  use arfima_kinds, only : dp
  use arfima_status, only : arfima_ok
  implicit none
  private

  integer, parameter, public :: long_memory_none = 0
  integer, parameter, public :: long_memory_fd = 1
  integer, parameter, public :: long_memory_fgn = 2
  integer, parameter, public :: long_memory_pla = 3

  type, public :: arfima_error
    integer :: code = arfima_ok
    character(len=:), allocatable :: message
  end type arfima_error

  type, public :: transfer_spec
    integer, allocatable :: r(:)
    integer, allocatable :: s(:)
    integer, allocatable :: b(:)
    real(dp), allocatable :: delta(:)
    real(dp), allocatable :: omega(:)
    real(dp), allocatable :: x(:,:)
  end type transfer_spec

  type, public :: arfima_spec
    integer :: p = 0
    integer :: q = 0
    integer :: pseas = 0
    integer :: qseas = 0
    integer :: dint = 0
    integer :: dseas = 0
    integer :: period = 0
    integer :: lmodel = long_memory_fd
    integer :: slmodel = long_memory_none
    logical :: estimate_mean = .true.
    logical :: use_regression = .false.
    logical :: use_transfer = .false.
    real(dp), allocatable :: xreg(:,:)
    type(transfer_spec) :: transfer
  end type arfima_spec

  type, public :: arfima_parameters
    real(dp), allocatable :: phi(:)
    real(dp), allocatable :: theta(:)
    real(dp), allocatable :: phiseas(:)
    real(dp), allocatable :: thetaseas(:)
    real(dp) :: dfrac = 0.0_dp
    real(dp) :: dfs = 0.0_dp
    real(dp) :: hurst = 0.5_dp
    real(dp) :: hurst_seasonal = 0.5_dp
    real(dp) :: alpha = 1.0_dp
    real(dp) :: alpha_seasonal = 1.0_dp
    real(dp) :: mean = 0.0_dp
    real(dp), allocatable :: beta(:)
    real(dp), allocatable :: delta(:)
    real(dp), allocatable :: omega(:)
  end type arfima_parameters

  type, public :: dl_result
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: sigma2_mle = huge(1.0_dp)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: innovation_variance(:)
    type(arfima_error) :: error
  end type dl_result

  type, public :: arfima_fit_result
    type(arfima_spec) :: spec
    type(arfima_parameters) :: parameters
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: sigma2 = huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: iterations = 0
    integer :: evaluations = 0
    logical :: converged = .false.
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: standard_error(:)
    type(arfima_error) :: error
  end type arfima_fit_result

  type, public :: arfima_forecast_result
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: standard_error(:)
    type(arfima_error) :: error
  end type arfima_forecast_result

  public :: set_error

  type, public :: arfima_mode_set
    type(arfima_fit_result), allocatable :: modes(:)
  end type arfima_mode_set

contains

  subroutine set_error(error, code, message)
    type(arfima_error), intent(out) :: error
    integer, intent(in) :: code
    character(len=*), intent(in) :: message
    error%code = code
    error%message = message
  end subroutine set_error

end module arfima_types
