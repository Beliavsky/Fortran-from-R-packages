! SPDX-License-Identifier: Artistic-2.0
module mts_types
   use mts_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: mts_success = 0
   integer, parameter, public :: mts_invalid_input = 1
   integer, parameter, public :: mts_singular = 2
   integer, parameter, public :: mts_no_convergence = 3

   type, public :: var_model
      integer :: order = 0
      integer :: n_series = 0
      logical :: include_mean = .true.
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: phi(:,:,:)
      real(dp), allocatable :: coef(:,:)
      real(dp), allocatable :: se_coef(:,:)
      real(dp), allocatable :: residuals(:,:)
      real(dp), allocatable :: sigma(:,:)
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp) :: hq = 0.0_dp
      integer :: n_parameters = 0
      integer :: status = mts_success
   end type var_model

   type, public :: varma_model
      integer :: p = 0
      integer :: q = 0
      integer :: n_series = 0
      logical :: include_mean = .true.
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: phi(:,:,:)
      real(dp), allocatable :: theta(:,:,:)
      real(dp), allocatable :: residuals(:,:)
      real(dp), allocatable :: sigma(:,:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = mts_success
   end type varma_model

   type, public :: varx_model
      integer :: p = 0
      integer :: exog_lags = 0
      integer :: n_series = 0
      integer :: n_exog = 0
      logical :: include_mean = .true.
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: phi(:,:,:)
      real(dp), allocatable :: beta(:,:,:)
      real(dp), allocatable :: residuals(:,:)
      real(dp), allocatable :: sigma(:,:)
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: status = mts_success
   end type varx_model

   type, public :: order_selection_result
      integer :: max_order = 0
      real(dp), allocatable :: aic(:)
      real(dp), allocatable :: bic(:)
      real(dp), allocatable :: hq(:)
      real(dp), allocatable :: m_stat(:)
      real(dp), allocatable :: p_value(:)
      integer :: aic_order = 0
      integer :: bic_order = 0
      integer :: hq_order = 0
      integer :: status = mts_success
   end type order_selection_result

   type, public :: diagnostic_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: degrees_freedom = 0
      integer :: status = mts_success
   end type diagnostic_result

   type, public :: dcc_model
      character(len=16) :: model_type = 'engle'
      character(len=16) :: distribution = 'normal'
      real(dp) :: alpha = 0.02_dp
      real(dp) :: beta = 0.95_dp
      real(dp) :: degrees_freedom = 8.0_dp
      real(dp), allocatable :: unconditional_corr(:,:)
      real(dp), allocatable :: correlations(:,:,:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = mts_success
   end type dcc_model

   type, public :: bekk_model
      logical :: include_mean = .true.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: c(:,:)
      real(dp), allocatable :: a(:,:)
      real(dp), allocatable :: b(:,:)
      real(dp), allocatable :: covariance(:,:,:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = mts_success
   end type bekk_model


   type, public :: vecm_model
      integer :: lag_order = 1
      integer :: rank = 0
      logical :: include_constant = .false.
      real(dp), allocatable :: alpha(:,:)
      real(dp), allocatable :: beta(:,:)
      real(dp), allocatable :: gamma(:,:,:)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: residuals(:,:)
      real(dp), allocatable :: sigma(:,:)
      real(dp), allocatable :: eigenvalues(:)
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: status = mts_success
   end type vecm_model

   type, public :: factor_result
      real(dp), allocatable :: loadings(:,:)
      real(dp), allocatable :: factors(:,:)
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: residuals(:,:)
      integer :: status = mts_success
   end type factor_result

   type, public :: bvar_result
      integer :: order = 0
      logical :: include_mean = .true.
      real(dp), allocatable :: posterior_mean(:,:)
      real(dp), allocatable :: posterior_covariance(:,:)
      real(dp), allocatable :: sigma(:,:)
      integer :: status = mts_success
   end type bvar_result

end module mts_types
