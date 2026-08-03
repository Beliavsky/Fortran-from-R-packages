! SPDX-License-Identifier: Artistic-2.0
module mts_rng
   use mts_kinds, only : dp
   use mts_linalg, only : cholesky_lower
   use mts_types, only : mts_success
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: set_random_seed, random_normal, random_chisq
   public :: random_multivariate_normal, random_multivariate_t
contains
   subroutine set_random_seed(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed+104729*i,huge(1)-1)
         if (put(i) <= 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine set_random_seed

   function random_normal() result(z)
      real(dp) :: z
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1,tiny(1.0_dp))
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function random_normal

   recursive function random_gamma(shape) result(x)
      real(dp), intent(in) :: shape
      real(dp) :: x, d, c, z, u
      if (shape <= 0.0_dp) then
         x = 0.0_dp
      else if (shape < 1.0_dp) then
         call random_number(u)
         x = random_gamma(shape+1.0_dp)*u**(1.0_dp/shape)
      else
         d = shape-1.0_dp/3.0_dp
         c = 1.0_dp/sqrt(9.0_dp*d)
         do
            z = random_normal()
            if (1.0_dp+c*z <= 0.0_dp) cycle
            x = d*(1.0_dp+c*z)**3
            call random_number(u)
            if (u < 1.0_dp-0.0331_dp*z**4) exit
            if (log(u) < 0.5_dp*z*z+d*(1.0_dp-x/d+log(x/d))) exit
         end do
      end if
   end function random_gamma

   function random_chisq(df) result(x)
      real(dp), intent(in) :: df
      real(dp) :: x
      x = 2.0_dp*random_gamma(0.5_dp*df)
   end function random_chisq

   subroutine random_multivariate_normal(mu,sigma,x,status)
      real(dp), intent(in) :: mu(:), sigma(:,:)
      real(dp), intent(out) :: x(size(mu))
      integer, intent(out), optional :: status
      real(dp), allocatable :: l(:,:)
      real(dp) :: z(size(mu))
      integer :: i, istat
      call cholesky_lower(sigma,l,istat,jitter=1.0e-12_dp)
      if (istat == mts_success) then
         do i = 1, size(mu)
            z(i) = random_normal()
         end do
         x = mu+matmul(l,z)
      else
         x = mu
      end if
      if (present(status)) status = istat
   end subroutine random_multivariate_normal

   subroutine random_multivariate_t(mu,sigma,df,x,variance_standardized,status)
      real(dp), intent(in) :: mu(:), sigma(:,:), df
      real(dp), intent(out) :: x(size(mu))
      logical, intent(in), optional :: variance_standardized
      integer, intent(out), optional :: status
      real(dp) :: z(size(mu)), scale
      integer :: istat
      logical :: standardized
      standardized = .false.
      if (present(variance_standardized)) standardized = variance_standardized
      call random_multivariate_normal(0.0_dp*mu,sigma,z,istat)
      if (istat == mts_success .and. df > 0.0_dp) then
         scale = sqrt(df/max(random_chisq(df),tiny(1.0_dp)))
         if (standardized .and. df > 2.0_dp) scale = scale*sqrt((df-2.0_dp)/df)
         x = mu+scale*z
      else
         x = mu
      end if
      if (present(status)) status = istat
   end subroutine random_multivariate_t
end module mts_rng
