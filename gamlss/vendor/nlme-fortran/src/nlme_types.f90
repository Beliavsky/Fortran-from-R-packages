! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_types
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS
  implicit none
  private

  integer, parameter, public :: COR_NONE = 0
  integer, parameter, public :: COR_AR1 = 1
  integer, parameter, public :: COR_CAR1 = 2
  integer, parameter, public :: COR_ARMA = 3
  integer, parameter, public :: COR_COMPOUND_SYMM = 4
  integer, parameter, public :: COR_EXPONENTIAL = 5
  integer, parameter, public :: COR_GAUSSIAN = 6
  integer, parameter, public :: COR_LINEAR = 7
  integer, parameter, public :: COR_RATIO = 8
  integer, parameter, public :: COR_SPHERICAL = 9
  integer, parameter, public :: COR_UNSTRUCTURED = 10

  integer, parameter, public :: VAR_CONSTANT = 0
  integer, parameter, public :: VAR_FIXED = 1
  integer, parameter, public :: VAR_IDENT = 2
  integer, parameter, public :: VAR_POWER = 3
  integer, parameter, public :: VAR_EXPONENTIAL = 4
  integer, parameter, public :: VAR_CONST_POWER = 5
  integer, parameter, public :: VAR_CONST_PROP = 6

  integer, parameter, public :: PD_IDENT = 0
  integer, parameter, public :: PD_DIAG = 1
  integer, parameter, public :: PD_LOG_CHOL = 2
  integer, parameter, public :: PD_COMPOUND_SYMM = 3

  type, public :: nlme_control
    integer :: max_iter = 250
    integer :: max_outer = 20
    real(dp) :: tolerance = 1.0e-7_dp
    real(dp) :: step = 0.15_dp
    real(dp) :: fd_step = 1.0e-5_dp
    logical :: reml = .true.
    logical :: optimize_covariance = .true.
    logical :: verbose = .false.
  end type nlme_control

  type, public :: correlation_spec
    integer :: kind = COR_NONE
    integer :: p = 0
    integer :: q = 0
    real(dp), allocatable :: par(:)
    logical :: nugget = .false.
    logical :: fixed = .false.
  end type correlation_spec

  type, public :: variance_spec
    integer :: kind = VAR_CONSTANT
    real(dp), allocatable :: par(:)
    logical :: fixed = .false.
  end type variance_spec

  type, public :: pd_spec
    integer :: kind = PD_LOG_CHOL
    integer :: dim = 0
    real(dp), allocatable :: par(:)
    logical :: fixed = .false.
  end type pd_spec

  type, public :: gls_result
    real(dp), allocatable :: beta(:)
    real(dp), allocatable :: beta_cov(:,:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: correlation_parameters(:)
    real(dp), allocatable :: variance_parameters(:)
    real(dp) :: sigma = 0.0_dp
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = NLME_SUCCESS
    logical :: converged = .false.
  end type gls_result

  type, public :: lme_result
    real(dp), allocatable :: beta(:)
    real(dp), allocatable :: beta_cov(:,:)
    real(dp), allocatable :: random_effects(:,:)
    integer, allocatable :: group_levels(:)
    real(dp), allocatable :: random_covariance(:,:)
    real(dp), allocatable :: fitted_marginal(:)
    real(dp), allocatable :: fitted_conditional(:)
    real(dp), allocatable :: residuals_marginal(:)
    real(dp), allocatable :: residuals_conditional(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: correlation_parameters(:)
    real(dp), allocatable :: variance_parameters(:)
    real(dp) :: sigma = 0.0_dp
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = NLME_SUCCESS
    logical :: converged = .false.
  end type lme_result

  type, public :: nonlinear_result
    real(dp), allocatable :: parameters(:)
    real(dp), allocatable :: parameter_cov(:,:)
    real(dp), allocatable :: random_effects(:,:)
    integer, allocatable :: group_levels(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residuals(:)
    real(dp) :: sigma = 0.0_dp
    real(dp) :: log_likelihood = -huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = NLME_SUCCESS
    logical :: converged = .false.
  end type nonlinear_result

  abstract interface
    subroutine nonlinear_model(theta, x, yhat, status)
      import dp
      real(dp), intent(in) :: theta(:)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: yhat(:)
      integer, intent(out) :: status
    end subroutine nonlinear_model
  end interface
  public :: nonlinear_model
end module nlme_types
