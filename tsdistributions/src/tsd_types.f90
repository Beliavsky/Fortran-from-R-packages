! SPDX-License-Identifier: GPL-2.0-only
module tsd_types
  use ghyp_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: tsd_success = 0
  integer, parameter, public :: tsd_invalid_argument = 1
  integer, parameter, public :: tsd_no_convergence = 2
  integer, parameter, public :: tsd_singular = 3
  integer, parameter, public :: tsd_numerical_failure = 4

  integer, parameter, public :: dist_norm = 1
  integer, parameter, public :: dist_std = 2
  integer, parameter, public :: dist_snorm = 3
  integer, parameter, public :: dist_sstd = 4
  integer, parameter, public :: dist_ged = 5
  integer, parameter, public :: dist_sged = 6
  integer, parameter, public :: dist_nig = 7
  integer, parameter, public :: dist_gh = 8
  integer, parameter, public :: dist_jsu = 9
  integer, parameter, public :: dist_ghst = 10

  type, public :: distribution_parameters
    real(dp) :: mu = 0.0_dp
    real(dp) :: sigma = 1.0_dp
    real(dp) :: skew = 1.0_dp
    real(dp) :: shape = 5.0_dp
    real(dp) :: lambda = -0.5_dp
  end type distribution_parameters

  type, public :: parameter_specification
    character(len=8) :: distribution = 'norm'
    type(distribution_parameters) :: parameters
    type(distribution_parameters) :: lower
    type(distribution_parameters) :: upper
    logical :: estimate_mu = .true.
    logical :: estimate_sigma = .true.
    logical :: estimate_skew = .false.
    logical :: estimate_shape = .false.
    logical :: estimate_lambda = .false.
  contains
    procedure :: number_estimated => specification_number_estimated
  end type parameter_specification

  type, public :: distribution_fit
    character(len=8) :: distribution = 'norm'
    type(distribution_parameters) :: parameters
    logical :: estimated(5) = [.true., .true., .false., .false., .false.]
    real(dp) :: negative_log_likelihood = huge(1.0_dp)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: hessian(:, :)
    real(dp), allocatable :: scores(:, :)
    real(dp), allocatable :: covariance_hessian(:, :)
    real(dp), allocatable :: covariance_opg(:, :)
    real(dp), allocatable :: covariance_qmle(:, :)
    integer :: nobs = 0
    integer :: degrees_of_freedom = 0
    integer :: iterations = 0
    integer :: status = tsd_invalid_argument
    character(len=200) :: message = ''
  end type distribution_fit

  type, public :: moment_summary
    real(dp) :: mean = 0.0_dp
    real(dp) :: standard_deviation = 1.0_dp
    real(dp) :: skewness = 0.0_dp
    real(dp) :: excess_kurtosis = 0.0_dp
    integer :: status = tsd_success
    character(len=160) :: message = ''
  end type moment_summary

  type, public :: profile_summary
    character(len=8) :: distribution = 'norm'
    integer, allocatable :: sizes(:)
    real(dp), allocatable :: rmse(:, :)
    real(dp), allocatable :: mae(:, :)
    real(dp), allocatable :: mape(:, :)
    real(dp) :: actual(7) = 0.0_dp
    integer :: successful_fits = 0
    integer :: attempted_fits = 0
    integer :: status = tsd_invalid_argument
    character(len=160) :: message = ''
  end type profile_summary

  type, public :: authorized_domain_result
    real(dp), allocatable :: skewness(:)
    real(dp), allocatable :: kurtosis(:)
    real(dp), allocatable :: skew_parameter(:)
    real(dp), allocatable :: shape_parameter(:)
    integer :: status = tsd_invalid_argument
    character(len=160) :: message = ''
  end type authorized_domain_result

  public :: distribution_id, valid_distribution, canonical_distribution_name

contains

  pure integer function specification_number_estimated(self) result(n)
    class(parameter_specification), intent(in) :: self
    n = count([self%estimate_mu, self%estimate_sigma, self%estimate_skew, &
               self%estimate_shape, self%estimate_lambda])
  end function specification_number_estimated

  pure integer function distribution_id(name) result(id)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: n
    n = canonical_distribution_name(name)
    select case (n)
    case ('norm');  id = dist_norm
    case ('std');   id = dist_std
    case ('snorm'); id = dist_snorm
    case ('sstd');  id = dist_sstd
    case ('ged');   id = dist_ged
    case ('sged');  id = dist_sged
    case ('nig');   id = dist_nig
    case ('gh');    id = dist_gh
    case ('jsu');   id = dist_jsu
    case ('ghst');  id = dist_ghst
    case default;   id = 0
    end select
  end function distribution_id

  pure logical function valid_distribution(name) result(ok)
    character(len=*), intent(in) :: name
    ok = distribution_id(name) > 0
  end function valid_distribution

  pure function canonical_distribution_name(name) result(value)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: value
    character(len=len_trim(name)) :: tmp
    integer :: i, code
    tmp = adjustl(name(1:len_trim(name)))
    do i = 1, len_trim(tmp)
      code = iachar(tmp(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) tmp(i:i) = achar(code + 32)
    end do
    value = trim(tmp)
  end function canonical_distribution_name

end module tsd_types
