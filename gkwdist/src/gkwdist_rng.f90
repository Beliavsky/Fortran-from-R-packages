! SPDX-License-Identifier: MIT
module gkwdist_rng
   use gkwdist_kinds, only : dp
   implicit none
   private
   public :: gamma_rng, beta_rng, seed_rng
contains
   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i=1,n
         put(i)=modulo(seed + 104729*i, huge(1)-1)
         if (put(i) == 0) put(i)=i
      end do
      call random_seed(put=put)
   end subroutine seed_rng

   function normal_rng() result(z)
      real(dp) :: z, u1, u2
      real(dp), parameter :: twopi=6.283185307179586476925286766559_dp
      call random_number(u1); call random_number(u2)
      u1=max(u1,tiny(1.0_dp))
      z=sqrt(-2.0_dp*log(u1))*cos(twopi*u2)
   end function normal_rng

   recursive function gamma_rng(shape) result(x)
      real(dp), intent(in) :: shape
      real(dp) :: x, d, c, z, v, u
      if (shape <= 0.0_dp) then
         x=0.0_dp; return
      end if
      if (shape < 1.0_dp) then
         call random_number(u)
         x=gamma_rng(shape+1.0_dp)*u**(1.0_dp/shape)
         return
      end if
      d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
      do
         do
            z=normal_rng(); v=1.0_dp+c*z
            if (v > 0.0_dp) exit
         end do
         v=v*v*v; call random_number(u)
         if (u < 1.0_dp-0.0331_dp*z**4) exit
         if (log(max(u,tiny(1.0_dp))) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
      end do
      x=d*v
   end function gamma_rng

   function beta_rng(a,b) result(x)
      real(dp), intent(in) :: a,b
      real(dp) :: x, ga, gb
      ga=gamma_rng(a); gb=gamma_rng(b)
      if (ga+gb > 0.0_dp) then
         x=ga/(ga+gb)
      else
         x=0.5_dp
      end if
   end function beta_rng
end module gkwdist_rng
