! Modern Fortran translation of R package skewunit.
! SPDX-License-Identifier: GPL-2.0-or-later
module skewunit_rng
   use skewunit_kinds, only : dp, nan_dp
   implicit none
   private

   public :: seed_skewunit_rng, rand_uniform, rand_normal, rand_gamma, rand_beta

contains

   subroutine seed_skewunit_rng(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)

      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729*i + 7919*i*i, huge(1)-1)
         if (put(i) == 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_skewunit_rng

   real(dp) function rand_uniform() result(u)
      call random_number(u)
      do while (u <= 0.0_dp .or. u >= 1.0_dp)
         call random_number(u)
      end do
   end function rand_uniform

   real(dp) function rand_normal() result(z)
      real(dp), save :: spare = 0.0_dp
      logical, save :: has_spare = .false.
      real(dp) :: u, v, s, fac

      if (has_spare) then
         z = spare
         has_spare = .false.
         return
      end if

      do
         u = 2.0_dp*rand_uniform()-1.0_dp
         v = 2.0_dp*rand_uniform()-1.0_dp
         s = u*u+v*v
         if (s > 0.0_dp .and. s < 1.0_dp) exit
      end do
      fac = sqrt(-2.0_dp*log(s)/s)
      z = u*fac
      spare = v*fac
      has_spare = .true.
   end function rand_normal

   recursive real(dp) function rand_gamma(shape, scale) result(x)
      real(dp), intent(in) :: shape, scale
      real(dp) :: d, c, z, v, u

      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         x = nan_dp()
         return
      end if

      if (shape < 1.0_dp) then
         x = rand_gamma(shape+1.0_dp,scale) &
             * rand_uniform()**(1.0_dp/shape)
         return
      end if

      d = shape-1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         do
            z = rand_normal()
            v = 1.0_dp+c*z
            if (v > 0.0_dp) exit
         end do
         v = v*v*v
         u = rand_uniform()
         if (u < 1.0_dp-0.0331_dp*z**4) exit
         if (log(u) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
      end do
      x = scale*d*v
   end function rand_gamma

   real(dp) function rand_beta(a, b) result(x)
      real(dp), intent(in) :: a, b
      real(dp) :: u, v

      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         x = nan_dp()
         return
      end if
      u = rand_gamma(a,1.0_dp)
      v = rand_gamma(b,1.0_dp)
      x = u/(u+v)
   end function rand_beta

end module skewunit_rng
