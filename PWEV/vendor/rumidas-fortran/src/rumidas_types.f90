! SPDX-License-Identifier: GPL-3.0-only
module rumidas_types
  use rumidas_kinds, only: dp
  use rumidas_status, only: RUMIDAS_INVALID_INPUT
  implicit none
  private

  integer, parameter, public :: RUMIDAS_GM = 1
  integer, parameter, public :: RUMIDAS_GM2M = 2
  integer, parameter, public :: RUMIDAS_GMX = 3
  integer, parameter, public :: RUMIDAS_DAGM = 4
  integer, parameter, public :: RUMIDAS_DAGM2M = 5
  integer, parameter, public :: RUMIDAS_DAGMX = 6

  integer, parameter, public :: RUMIDAS_MEM = 11
  integer, parameter, public :: RUMIDAS_MEM_MIDAS = 12
  integer, parameter, public :: RUMIDAS_MEM_X = 13
  integer, parameter, public :: RUMIDAS_MEM_MIDAS_X = 14

  integer, parameter, public :: RUMIDAS_NORMAL = 1
  integer, parameter, public :: RUMIDAS_STUDENT_T = 2
  integer, parameter, public :: RUMIDAS_BETA_LAG = 1
  integer, parameter, public :: RUMIDAS_ALMON_LAG = 2

  type, public :: garch_midas_spec
    integer :: model = RUMIDAS_GM
    integer :: distribution = RUMIDAS_NORMAL
    integer :: lag_function = RUMIDAS_BETA_LAG
    integer :: k1 = 1
    integer :: k2 = 0
    logical :: skew = .true.
  end type garch_midas_spec

  type, public :: mem_spec
    integer :: model = RUMIDAS_MEM
    integer :: k = 0
    logical :: skew = .true.
  end type mem_spec

  type, public :: rumidas_fit_control
    integer :: random_starts = 12
    integer :: max_iterations = 400
    integer :: random_seed = 12345
    character(len=16) :: method = 'bfgs'
    real(dp) :: gradient_tolerance = 1.0e-6_dp
    real(dp) :: relative_tolerance = 1.0e-9_dp
    logical :: compute_robust_covariance = .true.
  end type rumidas_fit_control

  type, public :: rumidas_fit_result
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: standard_errors(:)
    real(dp), allocatable :: robust_standard_errors(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: robust_covariance(:, :)
    real(dp), allocatable :: conditional(:)
    real(dp), allocatable :: long_run(:)
    real(dp), allocatable :: short_run(:)
    real(dp), allocatable :: loglik_obs(:)
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: iterations = 0
    integer :: function_count = 0
    integer :: status = RUMIDAS_INVALID_INPUT
    logical :: converged = .false.
    character(len=96) :: message = 'not run'
  end type rumidas_fit_result

end module rumidas_types
