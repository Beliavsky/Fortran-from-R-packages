! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_types
   use rmgarch_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: dist_gaussian = 1
   integer, parameter, public :: dist_student = 2
   integer, parameter, public :: dist_laplace = 3

   type, public :: univariate_garch_fit
      real(dp) :: mean = 0.0_dp
      real(dp) :: omega = 1.0e-6_dp
      real(dp) :: alpha = 0.05_dp
      real(dp) :: beta = 0.90_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: status = 1
      integer :: iterations = 0
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: standardized(:)
   end type univariate_garch_fit

   type, public :: dcc_spec
      integer :: p = 1
      integer :: q = 1
      integer :: g = 0
      integer :: distribution = dist_gaussian
      real(dp) :: shape = 8.0_dp
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: gamma(:)
   end type dcc_spec

   type, public :: dcc_fit_result
      type(dcc_spec) :: spec
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: status = 1
      integer :: iterations = 0
      character(len=128) :: message = 'not fitted'
      real(dp), allocatable :: qbar(:,:)
      real(dp), allocatable :: nbar(:,:)
      real(dp), allocatable :: q(:,:,:)
      real(dp), allocatable :: r(:,:,:)
      real(dp), allocatable :: loglikelihoods(:)
      real(dp), allocatable :: standardized_residuals(:,:)
   end type dcc_fit_result


   type, public :: multivariate_garch_fit
      type(univariate_garch_fit), allocatable :: margins(:)
      type(dcc_fit_result) :: dcc
      real(dp), allocatable :: sigma(:,:)
      real(dp), allocatable :: standardized(:,:)
      real(dp), allocatable :: covariance(:,:,:)
      integer :: status = 1
   end type multivariate_garch_fit

   type, public :: varx_fit_result
      integer :: order = 0
      logical :: include_intercept = .true.
      real(dp), allocatable :: coefficients(:,:)
      real(dp), allocatable :: residuals(:,:)
      real(dp), allocatable :: sigma(:,:)
      integer :: status = 1
   end type varx_fit_result

   type, public :: ica_result
      real(dp), allocatable :: mixing(:,:)
      real(dp), allocatable :: unmixing(:,:)
      real(dp), allocatable :: rotation(:,:)
      real(dp), allocatable :: sources(:,:)
      real(dp), allocatable :: center(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: whitening(:,:)
      real(dp), allocatable :: dewhitening(:,:)
      integer :: iterations = 0
      integer :: status = 1
   end type ica_result

   type, public :: gogarch_fit_result
      type(ica_result) :: ica
      type(univariate_garch_fit), allocatable :: components(:)
      real(dp), allocatable :: component_sigma(:,:)
      real(dp), allocatable :: standardized_components(:,:)
      real(dp), allocatable :: covariance(:,:,:)
      real(dp), allocatable :: correlation(:,:,:)
      integer :: status = 1
   end type gogarch_fit_result

   type, public :: fdcc_fit_result
      real(dp), allocatable :: group_alpha(:)
      real(dp), allocatable :: group_beta(:)
      integer, allocatable :: group_index(:)
      real(dp), allocatable :: a(:,:)
      real(dp), allocatable :: b(:,:)
      real(dp), allocatable :: c(:,:)
      real(dp), allocatable :: qbar(:,:)
      real(dp), allocatable :: q(:,:,:)
      real(dp), allocatable :: r(:,:,:)
      real(dp), allocatable :: loglikelihoods(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 1
      character(len=128) :: message = 'not fitted'
   end type fdcc_fit_result

   type, public :: copula_fit_result
      integer :: distribution = dist_gaussian
      logical :: time_varying = .false.
      real(dp) :: shape = 8.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: status = 1
      integer :: iterations = 0
      real(dp), allocatable :: correlation(:,:)
      type(dcc_fit_result) :: dcc
   end type copula_fit_result

   type, public :: copula_garch_fit_result
      type(univariate_garch_fit), allocatable :: margins(:)
      type(copula_fit_result) :: copula
      real(dp), allocatable :: sigma(:,:)
      real(dp), allocatable :: standardized(:,:)
      real(dp), allocatable :: uniforms(:,:)
      real(dp), allocatable :: scores(:,:)
      integer :: status = 1
   end type copula_garch_fit_result

   type, public :: rolling_gogarch_result
      integer :: window = 0
      integer :: refit_every = 1
      integer :: horizon = 1
      integer, allocatable :: origin(:)
      integer, allocatable :: fit_status(:)
      real(dp), allocatable :: covariance(:,:,:,:)
      real(dp), allocatable :: correlation(:,:,:,:)
   end type rolling_gogarch_result

   type, public :: rolling_dcc_result
      integer :: window = 0
      integer :: refit_every = 1
      integer :: horizon = 1
      integer, allocatable :: origin(:)
      integer, allocatable :: fit_status(:)
      real(dp), allocatable :: alpha(:,:)
      real(dp), allocatable :: beta(:,:)
      real(dp), allocatable :: gamma(:,:)
      real(dp), allocatable :: sigma(:,:,:)
      real(dp), allocatable :: covariance(:,:,:,:)
      real(dp), allocatable :: correlation(:,:,:,:)
   end type rolling_dcc_result

   type, public :: grid_distribution
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: density(:)
      real(dp), allocatable :: cdf(:)
      integer :: status = 1
   end type grid_distribution

contains

end module rmgarch_types
