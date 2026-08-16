! SPDX-License-Identifier: GPL-2.0-or-later
module fitdistrplus_types
  use fitdistrplus_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: fit_success = 0
  integer, parameter, public :: fit_invalid_argument = 1
  integer, parameter, public :: fit_no_convergence = 2
  integer, parameter, public :: fit_singular = 3
  integer, parameter, public :: fit_not_supported = 4
  integer, parameter, public :: fit_numerical_error = 5

  integer, parameter, public :: method_mle = 1
  integer, parameter, public :: method_mme = 2
  integer, parameter, public :: method_qme = 3
  integer, parameter, public :: method_mge = 4
  integer, parameter, public :: method_mse = 5

  integer, parameter, public :: gof_cvm = 1
  integer, parameter, public :: gof_ks = 2
  integer, parameter, public :: gof_ad = 3
  integer, parameter, public :: gof_adr = 4
  integer, parameter, public :: gof_adl = 5
  integer, parameter, public :: gof_ad2r = 6
  integer, parameter, public :: gof_ad2l = 7
  integer, parameter, public :: gof_ad2 = 8

  integer, parameter, public :: phi_kl = 1
  integer, parameter, public :: phi_j = 2
  integer, parameter, public :: phi_r = 3
  integer, parameter, public :: phi_h = 4
  integer, parameter, public :: phi_v = 5

  abstract interface
    function logpdf_callback(x, par) result(value)
      import dp
      real(dp), intent(in) :: x
      real(dp), intent(in) :: par(:)
      real(dp) :: value
    end function logpdf_callback

    function cdf_callback(x, par) result(value)
      import dp
      real(dp), intent(in) :: x
      real(dp), intent(in) :: par(:)
      real(dp) :: value
    end function cdf_callback

    function quantile_callback(prob, par) result(value)
      import dp
      real(dp), intent(in) :: prob
      real(dp), intent(in) :: par(:)
      real(dp) :: value
    end function quantile_callback

    function moment_callback(order, par) result(value)
      import dp
      integer, intent(in) :: order
      real(dp), intent(in) :: par(:)
      real(dp) :: value
    end function moment_callback

    function random_callback(par) result(value)
      import dp
      real(dp), intent(in) :: par(:)
      real(dp) :: value
    end function random_callback
  end interface

  type, public :: distribution_model
    character(len=32) :: name = "custom"
    integer :: npar = 0
    logical :: discrete = .false.
    character(len=24), allocatable :: parameter_names(:)
    real(dp), allocatable :: default_lower(:)
    real(dp), allocatable :: default_upper(:)
    procedure(logpdf_callback), pointer, nopass :: logpdf => null()
    procedure(cdf_callback), pointer, nopass :: cdf => null()
    procedure(quantile_callback), pointer, nopass :: quantile => null()
    procedure(moment_callback), pointer, nopass :: raw_moment => null()
    procedure(random_callback), pointer, nopass :: random_value => null()
  end type distribution_model

  type, public :: fit_control
    integer :: max_iterations = 2000
    real(dp) :: tolerance = 1.0e-8_dp
    real(dp) :: simplex_scale = 0.10_dp
    logical :: calculate_vcov = .true.
  end type fit_control

  type, public :: fit_result
    character(len=32) :: distribution = ""
    character(len=16) :: method = ""
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: standard_error(:)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: nobs = 0
    integer :: iterations = 0
    integer :: convergence = fit_invalid_argument
    character(len=160) :: message = "not fitted"
  end type fit_result

  type, public :: censored_sample
    real(dp), allocatable :: left(:)
    real(dp), allocatable :: right(:)
  end type censored_sample

  type, public :: descriptive_result
    integer :: n = 0
    real(dp) :: minimum = 0.0_dp
    real(dp) :: maximum = 0.0_dp
    real(dp) :: median = 0.0_dp
    real(dp) :: mean = 0.0_dp
    real(dp) :: standard_deviation = 0.0_dp
    real(dp) :: skewness = 0.0_dp
    real(dp) :: kurtosis = 0.0_dp
    logical :: unbiased = .true.
    integer :: status = fit_invalid_argument
  end type descriptive_result

  type, public :: gof_result
    real(dp) :: ks = huge(1.0_dp)
    real(dp) :: cvm = huge(1.0_dp)
    real(dp) :: ad = huge(1.0_dp)
    real(dp) :: chi_square = huge(1.0_dp)
    real(dp) :: chi_square_pvalue = -1.0_dp
    integer :: chi_square_df = 0
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: status = fit_invalid_argument
  end type gof_result

  type, public :: bootstrap_result
    real(dp), allocatable :: estimates(:, :)
    integer, allocatable :: convergence(:)
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: standard_deviation(:)
    integer :: nsim = 0
    integer :: successful = 0
    integer :: status = fit_invalid_argument
  end type bootstrap_result

  type, public :: npmle_result
    real(dp), allocatable :: interval_left(:)
    real(dp), allocatable :: interval_right(:)
    real(dp), allocatable :: probability(:)
    integer :: iterations = 0
    integer :: status = fit_invalid_argument
  end type npmle_result

end module fitdistrplus_types
