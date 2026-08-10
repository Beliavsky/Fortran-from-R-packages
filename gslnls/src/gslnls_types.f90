! SPDX-License-Identifier: LGPL-3.0-only
module gslnls_types
  use gslnls_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: NLS_SUCCESS = 0
  integer, parameter, public :: NLS_MAXITER = 1
  integer, parameter, public :: NLS_NO_PROGRESS = 2
  integer, parameter, public :: NLS_BAD_FUNCTION = 3
  integer, parameter, public :: NLS_BAD_INPUT = 4
  integer, parameter, public :: NLS_SINGULAR = 5

  integer, parameter, public :: NLS_LM = 1
  integer, parameter, public :: NLS_LMACCEL = 2
  integer, parameter, public :: NLS_DOGLEG = 3
  integer, parameter, public :: NLS_DDOGLEG = 4
  integer, parameter, public :: NLS_SUBSPACE2D = 5
  integer, parameter, public :: NLS_CGST = 6

  integer, parameter, public :: NLS_SCALE_MORE = 1
  integer, parameter, public :: NLS_SCALE_LEVENBERG = 2
  integer, parameter, public :: NLS_SCALE_MARQUARDT = 3

  integer, parameter, public :: NLS_SOLVER_QR = 1
  integer, parameter, public :: NLS_SOLVER_CHOLESKY = 2
  integer, parameter, public :: NLS_SOLVER_SVD = 3

  integer, parameter, public :: NLS_FD_FORWARD = 1
  integer, parameter, public :: NLS_FD_CENTER = 2

  integer, parameter, public :: LOSS_DEFAULT = 0
  integer, parameter, public :: LOSS_HUBER = 1
  integer, parameter, public :: LOSS_BARRON = 2
  integer, parameter, public :: LOSS_BISQUARE = 3
  integer, parameter, public :: LOSS_WELSH = 4
  integer, parameter, public :: LOSS_OPTIMAL = 5
  integer, parameter, public :: LOSS_HAMPEL = 6
  integer, parameter, public :: LOSS_GGW = 7
  integer, parameter, public :: LOSS_LQQ = 8

  type, public :: nls_control
    integer :: maxiter = 100
    integer :: algorithm = NLS_LM
    integer :: scale = NLS_SCALE_MORE
    integer :: solver = NLS_SOLVER_QR
    integer :: fdtype = NLS_FD_FORWARD
    real(dp) :: factor_up = 2.0_dp
    real(dp) :: factor_down = 3.0_dp
    real(dp) :: avmax = 0.75_dp
    real(dp) :: h_df = sqrt(epsilon(1.0_dp))
    real(dp) :: h_fvv = 0.02_dp
    real(dp) :: xtol = sqrt(epsilon(1.0_dp))
    real(dp) :: ftol = sqrt(epsilon(1.0_dp))
    real(dp) :: gtol = sqrt(epsilon(1.0_dp))
    integer :: mstart_n = 30
    integer :: mstart_p = 5
    integer :: mstart_q = 3
    real(dp) :: mstart_r = 4.0_dp
    integer :: mstart_s = 2
    real(dp) :: mstart_tol = 0.25_dp
    integer :: mstart_maxiter = 10
    integer :: mstart_maxstart = 250
    integer :: mstart_minsp = 1
    integer :: irls_maxiter = 50
    real(dp) :: irls_xtol = epsilon(1.0_dp)**0.25_dp
    logical :: store_trace = .false.
  end type nls_control

  type, public :: nls_loss
    integer :: kind = LOSS_DEFAULT
    real(dp) :: cc(3) = 0.0_dp
    integer :: ncc = 0
  end type nls_loss

  type, public :: nls_result
    real(dp), allocatable :: par(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residual(:)
    real(dp), allocatable :: weighted_residual(:)
    real(dp), allocatable :: jacobian(:,:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: irls_weights(:)
    real(dp), allocatable :: irls_psi(:)
    real(dp), allocatable :: irls_dpsi(:)
    real(dp), allocatable :: par_trace(:,:)
    real(dp), allocatable :: ssr_trace(:)
    real(dp) :: ssr = huge(1.0_dp)
    real(dp) :: sigma = huge(1.0_dp)
    real(dp) :: irls_sigma = 0.0_dp
    real(dp) :: gradient_inf = huge(1.0_dp)
    integer :: iterations = 0
    integer :: evaluations = 0
    integer :: jacobian_evaluations = 0
    integer :: fvv_evaluations = 0
    integer :: irls_iterations = 0
    integer :: rank = 0
    integer :: status = NLS_BAD_INPUT
    integer :: info = 0
    logical :: converged = .false.
  end type nls_result

  type, public :: multistart_result
    type(nls_result) :: fit
    integer :: starts_evaluated = 0
    integer :: local_searches = 0
    integer :: major_iterations = 0
    real(dp), allocatable :: best_start(:)
  end type multistart_result

  abstract interface
    subroutine nls_model(par, yhat, ierr)
      import dp
      real(dp), intent(in) :: par(:)
      real(dp), intent(out) :: yhat(:)
      integer, intent(out) :: ierr
    end subroutine nls_model

    subroutine nls_jacobian(par, jac, ierr)
      import dp
      real(dp), intent(in) :: par(:)
      real(dp), intent(out) :: jac(:,:)
      integer, intent(out) :: ierr
    end subroutine nls_jacobian

    subroutine nls_fvv(par, v, fvv, ierr)
      import dp
      real(dp), intent(in) :: par(:), v(:)
      real(dp), intent(out) :: fvv(:)
      integer, intent(out) :: ierr
    end subroutine nls_fvv

    subroutine nls_jacobian_operator(par, transpose_j, u, v, ierr)
      import dp
      real(dp), intent(in) :: par(:), u(:)
      logical, intent(in) :: transpose_j
      real(dp), intent(out) :: v(:)
      integer, intent(out) :: ierr
    end subroutine nls_jacobian_operator
  end interface

  public :: nls_model, nls_jacobian, nls_fvv, nls_jacobian_operator
  public :: default_loss

contains

  pure function default_loss(kind) result(loss)
    integer, intent(in), optional :: kind
    type(nls_loss) :: loss
    integer :: k

    k = LOSS_DEFAULT
    if (present(kind)) k = kind
    loss%kind = k
    select case (k)
    case (LOSS_DEFAULT)
      loss%ncc = 0
    case (LOSS_HUBER)
      loss%ncc = 1; loss%cc(1) = 1.345_dp
    case (LOSS_BARRON)
      loss%ncc = 2; loss%cc(1:2) = [1.0_dp, 1.345_dp]
    case (LOSS_BISQUARE)
      loss%ncc = 1; loss%cc(1) = 4.685061_dp
    case (LOSS_WELSH)
      loss%ncc = 1; loss%cc(1) = 2.11_dp
    case (LOSS_OPTIMAL)
      loss%ncc = 1; loss%cc(1) = 1.060158_dp
    case (LOSS_HAMPEL)
      loss%ncc = 1; loss%cc(1) = 0.9016085_dp
    case (LOSS_GGW)
      loss%ncc = 3; loss%cc = [1.387_dp, 1.5_dp, 1.063_dp]
    case (LOSS_LQQ)
      loss%ncc = 3; loss%cc = [1.473_dp, 0.982_dp, 1.5_dp]
    case default
      loss%kind = LOSS_DEFAULT
      loss%ncc = 0
    end select
  end function default_loss

end module gslnls_types
