! SPDX-License-Identifier: GPL-3.0-only
module ufrisk_types
   use kind_mod, only : dp
   use rugarch_mod, only : rugarch_fit_t
   use smoots_types, only : arma_model, smooth_result
   use fracdiff_types, only : fracdiff_model
   implicit none
   private

   integer, parameter, public :: ufrisk_model_sgarch = 1
   integer, parameter, public :: ufrisk_model_lgarch = 2
   integer, parameter, public :: ufrisk_model_egarch = 3
   integer, parameter, public :: ufrisk_model_aparch = 4
   integer, parameter, public :: ufrisk_model_figarch = 5
   integer, parameter, public :: ufrisk_model_filgarch = 6

   integer, parameter, public :: ufrisk_distribution_normal = 1
   integer, parameter, public :: ufrisk_distribution_student = 2
   integer, parameter, public :: ufrisk_smooth_none = 0
   integer, parameter, public :: ufrisk_smooth_lpr = 1

   integer, parameter, public :: ufrisk_ok = 0
   integer, parameter, public :: ufrisk_invalid_input = 1
   integer, parameter, public :: ufrisk_smoothing_failed = 2
   integer, parameter, public :: ufrisk_model_fit_failed = 3
   integer, parameter, public :: ufrisk_numerical_failure = 4
   integer, parameter, public :: ufrisk_no_violations = 5

   type, public :: ufrisk_options
      integer :: model = ufrisk_model_sgarch
      integer :: distribution = ufrisk_distribution_student
      integer :: smooth = ufrisk_smooth_none
      integer :: arch_order = 1
      integer :: garch_order = 1
      integer :: n_out = 250
      integer :: smoothing_order = 3
      integer :: smoothing_mu = 1
      integer :: smoothing_iterations = 40
      integer :: fractional_p_min = 0
      integer :: fractional_p_max = 0
      integer :: fractional_q_min = 0
      integer :: fractional_q_max = 0
      integer :: fractional_terms = 100
      integer :: truncation_lag = 1000
      integer :: max_fit_iterations = 250
      integer :: log_filter_lag = 50
      real(dp) :: var_confidence = 0.99_dp
      real(dp) :: es_confidence = 0.975_dp
      real(dp) :: smoothing_start = 0.15_dp
      real(dp) :: fit_tolerance = 1.0e-5_dp
      character(len=8) :: smoothing_algorithm = 'A'
   end type ufrisk_options

   type, public :: long_memory_smooth_result
      integer :: status = ufrisk_ok
      integer :: iterations = 0
      integer :: ar_order = 0
      integer :: ma_order = 0
      real(dp) :: bandwidth = 0.0_dp
      real(dp) :: d = 0.0_dp
      real(dp) :: cf0 = 0.0_dp
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: bandwidth_path(:)
      type(fracdiff_model) :: fractional_fit
   end type long_memory_smooth_result

   type, public :: ufrisk_result
      integer :: status = ufrisk_ok
      integer :: model = ufrisk_model_sgarch
      integer :: distribution = ufrisk_distribution_student
      integer :: smooth = ufrisk_smooth_none
      integer :: arch_order = 1
      integer :: garch_order = 1
      real(dp) :: mean_return = 0.0_dp
      real(dp) :: degrees_freedom = 0.0_dp
      real(dp) :: var_tail_probability = 0.01_dp
      real(dp) :: es_tail_probability = 0.025_dp
      real(dp) :: scale_forecast = 1.0_dp
      character(len=24) :: model_name = 'sGARCH'
      character(len=160) :: message = ''
      real(dp), allocatable :: returns_in(:)
      real(dp), allocatable :: returns_out(:)
      real(dp), allocatable :: centered_in(:)
      real(dp), allocatable :: centered_out(:)
      real(dp), allocatable :: sigma_in(:)
      real(dp), allocatable :: sigma_forecast(:)
      real(dp), allocatable :: scale(:)
      real(dp), allocatable :: var_es_level(:)
      real(dp), allocatable :: var_var_level(:)
      real(dp), allocatable :: expected_shortfall(:)
      type(rugarch_fit_t) :: rugarch_fit
      type(arma_model) :: arma_fit
      type(fracdiff_model) :: fracdiff_fit
      type(smooth_result) :: short_memory_smooth
      type(long_memory_smooth_result) :: long_memory_smooth
   end type ufrisk_result

   type, public :: ufrisk_traffic_result
      integer :: status = ufrisk_ok
      integer :: violations_var = 0
      integer :: violations_es_var = 0
      integer :: violations_es = 0
      real(dp) :: p_var = 1.0_dp
      real(dp) :: p_es_var = 1.0_dp
      real(dp) :: p_es = 1.0_dp
      real(dp) :: breach_sum = 0.0_dp
      real(dp) :: weighted_absolute_deviation = 0.0_dp
   end type ufrisk_traffic_result

   type, public :: ufrisk_coverage_result
      integer :: status = ufrisk_ok
      real(dp) :: tail_probability = 0.0_dp
      real(dp) :: lr_unconditional = 0.0_dp
      real(dp) :: lr_independence = 0.0_dp
      real(dp) :: lr_conditional = 0.0_dp
      real(dp) :: p_unconditional = 1.0_dp
      real(dp) :: p_independence = 1.0_dp
      real(dp) :: p_conditional = 1.0_dp
      integer :: n00 = 0
      integer :: n01 = 0
      integer :: n10 = 0
      integer :: n11 = 0
   end type ufrisk_coverage_result

   type, public :: ufrisk_loss_result
      integer :: status = ufrisk_ok
      real(dp) :: regulatory = 0.0_dp
      real(dp) :: firm = 0.0_dp
      real(dp) :: abad = 0.0_dp
      real(dp) :: feng = 0.0_dp
   end type ufrisk_loss_result

   public :: ufrisk_model_label, ufrisk_status_message
contains
   pure function ufrisk_model_label(model, semiparametric) result(label)
      integer, intent(in) :: model
      logical, intent(in), optional :: semiparametric
      character(len=24) :: label
      logical :: semi
      semi = .false.
      if (present(semiparametric)) semi = semiparametric
      select case (model)
      case (ufrisk_model_sgarch)
         label = 'sGARCH'
      case (ufrisk_model_lgarch)
         label = 'lGARCH'
      case (ufrisk_model_egarch)
         label = 'eGARCH'
      case (ufrisk_model_aparch)
         label = 'apARCH'
      case (ufrisk_model_figarch)
         label = 'fiGARCH'
      case (ufrisk_model_filgarch)
         label = 'filGARCH'
      case default
         label = 'unknown'
      end select
      if (semi) label = 'Semi-' // trim(label)
   end function ufrisk_model_label

   pure function ufrisk_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=80) :: message
      select case (status)
      case (ufrisk_ok)
         message = 'success'
      case (ufrisk_invalid_input)
         message = 'invalid input'
      case (ufrisk_smoothing_failed)
         message = 'nonparametric scale estimation failed'
      case (ufrisk_model_fit_failed)
         message = 'parametric model estimation failed'
      case (ufrisk_numerical_failure)
         message = 'numerical calculation failed'
      case (ufrisk_no_violations)
         message = 'no VaR violations; coverage tests are not applicable'
      case default
         message = 'unknown status'
      end select
   end function ufrisk_status_message
end module ufrisk_types
