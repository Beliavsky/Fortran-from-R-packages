! SPDX-License-Identifier: GPL-2.0-or-later
module fints_types
   use fints_kinds, only : dp
   use fints_status, only : fints_ok
   implicit none
   private

   type, public :: test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: degrees_freedom = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: lag = 0
      integer :: n_observations = 0
      integer :: status = fints_ok
      character(len=64) :: method = ''
   end type test_result

   type, public :: summary_result
      real(dp) :: start = 1.0_dp
      integer :: size = 0
      real(dp) :: mean = 0.0_dp
      real(dp) :: standard_deviation = 0.0_dp
      real(dp) :: skewness = 0.0_dp
      real(dp) :: excess_kurtosis = 0.0_dp
      real(dp) :: minimum = 0.0_dp
      real(dp) :: maximum = 0.0_dp
      integer :: status = fints_ok
   end type summary_result

   type, public :: acf_result
      real(dp), allocatable :: lag(:)
      real(dp), allocatable :: value(:)
      integer :: n_used = 0
      integer :: status = fints_ok
      character(len=16) :: acf_type = 'correlation'
   end type acf_result

   type, public :: cross_acf_result
      real(dp), allocatable :: lag(:)
      real(dp), allocatable :: value(:,:,:)
      integer :: n_used = 0
      integer :: status = fints_ok
      character(len=16) :: acf_type = 'correlation'
   end type cross_acf_result

   type, public :: apca_result
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: factors(:,:)
      real(dp), allocatable :: loadings(:,:)
      real(dp), allocatable :: r_squared(:)
      integer :: status = fints_ok
   end type apca_result

   type, public :: arma_acf_result
      complex(dp), allocatable :: roots(:)
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: damping(:)
      real(dp), allocatable :: period(:)
      logical :: stationary = .true.
      logical :: partial = .false.
      integer :: status = fints_ok
   end type arma_acf_result

   type, public :: yearmon_result
      integer, allocatable :: year(:)
      integer, allocatable :: month(:)
      logical :: converted = .false.
      integer :: duplicate_count = 0
      integer :: status = fints_ok
   end type yearmon_result

   type, public :: arima_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
      real(dp), allocatable :: seasonal_ar(:)
      real(dp), allocatable :: seasonal_ma(:)
      real(dp), allocatable :: regression(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: fitted(:)
      real(dp) :: intercept = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: r_squared = 0.0_dp
      integer :: n_used = 0
      integer :: iterations = 0
      logical :: converged = .false.
      integer :: status = fints_ok
      character(len=16) :: method = 'CSS-ML'
      type(test_result) :: box_test
   end type arima_result
end module fints_types
