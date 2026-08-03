! SPDX-License-Identifier: GPL-3.0-only
module fitheavytail_types
   use fitheavytail_kinds, only: dp
   use fitheavytail_status, only: ht_success
   implicit none
   private

   public :: clear_fit

   type, public :: heavy_tail_fit
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: gamma(:)
      real(dp), allocatable :: scatter(:,:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: loadings(:,:)
      real(dp), allocatable :: psi(:)
      real(dp), allocatable :: latent_weights(:)
      real(dp) :: nu = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: cpu_time = 0.0_dp
      logical :: converged = .false.
      integer :: num_iterations = 0
      integer :: status = ht_success
      character(len=160) :: message = ''
   end type heavy_tail_fit

contains

   subroutine clear_fit(result)
      type(heavy_tail_fit), intent(inout) :: result
      if (allocated(result%mu)) deallocate(result%mu)
      if (allocated(result%gamma)) deallocate(result%gamma)
      if (allocated(result%scatter)) deallocate(result%scatter)
      if (allocated(result%covariance)) deallocate(result%covariance)
      if (allocated(result%mean)) deallocate(result%mean)
      if (allocated(result%loadings)) deallocate(result%loadings)
      if (allocated(result%psi)) deallocate(result%psi)
      if (allocated(result%latent_weights)) deallocate(result%latent_weights)
      result%nu = 0.0_dp
      result%log_likelihood = 0.0_dp
      result%cpu_time = 0.0_dp
      result%converged = .false.
      result%num_iterations = 0
      result%status = ht_success
      result%message = ''
   end subroutine clear_fit

end module fitheavytail_types
