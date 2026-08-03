! SPDX-License-Identifier: GPL-3.0-only
module highorder_types
   use fitheavytail_kinds, only: dp
   implicit none
   private

   integer, parameter, public :: hop_success = 0
   integer, parameter, public :: hop_invalid_argument = 1
   integer, parameter, public :: hop_dimension_mismatch = 2
   integer, parameter, public :: hop_numerical_error = 3
   integer, parameter, public :: hop_not_converged = 4
   integer, parameter, public :: hop_fit_error = 5

   type, public :: sample_moments
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: centered(:,:)
      real(dp), allocatable :: coskewness(:,:,:)
      real(dp), allocatable :: cokurtosis(:,:,:,:)
      real(dp) :: third_scale = 1.0_dp
      real(dp) :: fourth_scale = 1.0_dp
      integer :: nobs = 0
      integer :: nassets = 0
      logical :: magnitude_adjusted = .false.
      integer :: status = hop_success
      character(len=200) :: message = ''
   end type sample_moments

   type, public :: skew_t_parameters
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: gamma(:)
      real(dp), allocatable :: scatter(:,:)
      real(dp), allocatable :: chol_scatter(:,:)
      real(dp) :: nu = 0.0_dp
      real(dp) :: a11 = 0.0_dp
      real(dp) :: a21 = 0.0_dp
      real(dp) :: a22 = 0.0_dp
      real(dp) :: a31 = 0.0_dp
      real(dp) :: a32 = 0.0_dp
      real(dp) :: a41 = 0.0_dp
      real(dp) :: a42 = 0.0_dp
      real(dp) :: a43 = 0.0_dp
      integer :: num_iterations = 0
      logical :: converged = .false.
      integer :: status = hop_success
      character(len=200) :: message = ''
   end type skew_t_parameters

   type, public :: portfolio_result
      real(dp), allocatable :: w(:)
      real(dp), allocatable :: objective_history(:)
      real(dp) :: moments(4) = 0.0_dp
      real(dp) :: delta = 0.0_dp
      real(dp) :: improvement(4) = 0.0_dp
      integer :: iterations = 0
      logical :: converged = .false.
      integer :: status = hop_success
      character(len=200) :: message = ''
   end type portfolio_result

   public :: clear_sample_moments, clear_skew_t_parameters, clear_portfolio_result

contains

   subroutine clear_sample_moments(x)
      type(sample_moments), intent(inout) :: x
      if (allocated(x%mu)) deallocate(x%mu)
      if (allocated(x%covariance)) deallocate(x%covariance)
      if (allocated(x%centered)) deallocate(x%centered)
      if (allocated(x%coskewness)) deallocate(x%coskewness)
      if (allocated(x%cokurtosis)) deallocate(x%cokurtosis)
      x%third_scale = 1.0_dp
      x%fourth_scale = 1.0_dp
      x%nobs = 0
      x%nassets = 0
      x%magnitude_adjusted = .false.
      x%status = hop_success
      x%message = ''
   end subroutine clear_sample_moments

   subroutine clear_skew_t_parameters(x)
      type(skew_t_parameters), intent(inout) :: x
      if (allocated(x%mu)) deallocate(x%mu)
      if (allocated(x%gamma)) deallocate(x%gamma)
      if (allocated(x%scatter)) deallocate(x%scatter)
      if (allocated(x%chol_scatter)) deallocate(x%chol_scatter)
      x%nu = 0.0_dp
      x%a11 = 0.0_dp
      x%a21 = 0.0_dp
      x%a22 = 0.0_dp
      x%a31 = 0.0_dp
      x%a32 = 0.0_dp
      x%a41 = 0.0_dp
      x%a42 = 0.0_dp
      x%a43 = 0.0_dp
      x%num_iterations = 0
      x%converged = .false.
      x%status = hop_success
      x%message = ''
   end subroutine clear_skew_t_parameters

   subroutine clear_portfolio_result(x)
      type(portfolio_result), intent(inout) :: x
      if (allocated(x%w)) deallocate(x%w)
      if (allocated(x%objective_history)) deallocate(x%objective_history)
      x%moments = 0.0_dp
      x%delta = 0.0_dp
      x%improvement = 0.0_dp
      x%iterations = 0
      x%converged = .false.
      x%status = hop_success
      x%message = ''
   end subroutine clear_portfolio_result

end module highorder_types
