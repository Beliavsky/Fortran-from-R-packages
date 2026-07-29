! SPDX-License-Identifier: GPL-3.0-only
module smoots_types
   use smoots_kinds, only : dp
   use smoots_status, only : sm_ok
   implicit none
   private

   integer, parameter, public :: sm_method_lpr = 1
   integer, parameter, public :: sm_method_kernel = 2
   integer, parameter, public :: sm_cf_np = 1
   integer, parameter, public :: sm_cf_arma = 2
   integer, parameter, public :: sm_cf_ar = 3
   integer, parameter, public :: sm_cf_ma = 4
   integer, parameter, public :: sm_infl_opt = 1
   integer, parameter, public :: sm_infl_naive = 2
   integer, parameter, public :: sm_infl_variance = 3
   integer, parameter, public :: sm_trend_linear = 1
   integer, parameter, public :: sm_trend_constant = 2

   type, public :: arma_model
      integer :: p = 0
      integer :: q = 0
      logical :: include_mean = .false.
      integer :: status = sm_ok
      integer :: iterations = 0
      real(dp) :: mean_value = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: fitted(:)
   end type arma_model

   type, public :: smooth_result
      integer :: status = sm_ok
      integer :: n = 0
      integer :: v = 0
      integer :: p = 1
      integer :: mu = 1
      integer :: bb = 1
      integer :: method = sm_method_lpr
      integer :: cf_method = sm_cf_np
      integer :: inflation = sm_infl_opt
      integer :: pilot_p = 1
      integer :: niterations = 0
      integer :: lag_window = 0
      integer :: ar_order = 0
      integer :: ma_order = 0
      logical :: enlarged_variance_bandwidth = .true.
      real(dp) :: b0 = 0.0_dp
      real(dp) :: b_start = 0.15_dp
      real(dp) :: pilot_b_start = 0.15_dp
      real(dp) :: cf0 = 0.0_dp
      real(dp) :: curvature_integral = 0.0_dp
      real(dp) :: boundary_cut = 0.05_dp
      real(dp), allocatable :: original(:)
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: weights(:,:)
      real(dp), allocatable :: bandwidth_steps(:)
   end type smooth_result

   type, public :: confidence_result
      integer :: status = sm_ok
      integer :: derivative_order = 0
      real(dp) :: confidence_level = 0.95_dp
      real(dp) :: unbiased_bandwidth = 0.0_dp
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp), allocatable :: parametric(:)
   end type confidence_result

   type, public :: forecast_result
      integer :: status = sm_ok
      integer :: horizon = 0
      integer :: simulations = 0
      real(dp) :: confidence_level = 0.95_dp
      real(dp), allocatable :: point(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp), allocatable :: errors(:,:)
      type(arma_model) :: model
   end type forecast_result

   type, public :: rolling_result
      integer :: status = sm_ok
      integer :: horizon = 0
      real(dp) :: mase = 0.0_dp
      real(dp) :: rmsse = 0.0_dp
      real(dp), allocatable :: observed(:)
      real(dp), allocatable :: trend_forecast(:)
      real(dp), allocatable :: residual_forecast(:)
      real(dp), allocatable :: point(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      logical, allocatable :: breach(:)
      real(dp), allocatable :: breach_value(:)
      type(smooth_result) :: smooth_model
      type(arma_model) :: arma
   end type rolling_result
end module smoots_types
