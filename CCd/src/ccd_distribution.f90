module ccd_distribution
   use ccd_kinds, only : dp, i8
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: dcc, pcc, qcc
contains
   elemental function dcc(y, mu, lambda, logged) result(ans)
      real(dp), intent(in) :: y
      real(dp), intent(in), optional :: mu
      real(dp), intent(in) :: lambda
      logical, intent(in), optional :: logged
      real(dp) :: ans, loc, logden
      logical :: want_log

      loc = 0.0_dp
      if (present(mu)) loc = mu
      want_log = .false.
      if (present(logged)) want_log = logged

      if (lambda <= 0.0_dp) then
         ans = ieee_nan()
         return
      end if

      logden = log(tanh(lambda*pi)) - log(pi) + log(lambda) &
         - log(lambda*lambda + (y-loc)*(y-loc))
      if (want_log) then
         ans = logden
      else
         ans = exp(logden)
      end if
   end function dcc

   elemental function pcc(y, mu, lambda) result(ans)
      integer(i8), intent(in) :: y
      real(dp), intent(in), optional :: mu
      real(dp), intent(in) :: lambda
      real(dp) :: ans, loc, s
      integer(i8) :: k

      loc = 0.0_dp
      if (present(mu)) loc = mu
      if (lambda <= 0.0_dp) then
         ans = ieee_nan()
         return
      end if

      if (y == 0_i8) then
         ans = 0.5_dp + 0.5_dp*dcc(0.0_dp, loc, lambda)
      else if (y < 0_i8) then
         s = 0.0_dp
         ! This reproduces the upstream R sequence (y+1):(abs(y)-1).
         do k = y + 1_i8, abs(y) - 1_i8
            s = s + dcc(real(k, dp), loc, lambda)
         end do
         ans = 0.5_dp - 0.5_dp*s
      else
         s = 0.0_dp
         do k = 1_i8, y
            s = s + dcc(real(k, dp), loc, lambda)
         end do
         ans = 0.5_dp + 0.5_dp*dcc(0.0_dp, loc, lambda) + s
      end if
   end function pcc

   elemental function qcc(p, mu, lambda) result(q)
      real(dp), intent(in) :: p, mu, lambda
      integer(i8) :: q
      integer(i8) :: lo, hi, mid
      real(dp) :: c

      if (p <= 0.0_dp) then
         q = -huge(0_i8)
         return
      else if (p >= 1.0_dp) then
         q = huge(0_i8)
         return
      else if (lambda <= 0.0_dp) then
         q = 0_i8
         return
      end if

      lo = -1_i8
      hi = 1_i8
      do while (pcc(lo, mu, lambda) >= p)
         if (lo < -1000000000_i8) exit
         lo = 2_i8*lo
      end do
      do while (pcc(hi, mu, lambda) < p)
         if (hi > 1000000000_i8) exit
         hi = 2_i8*hi
      end do

      do while (hi - lo > 1_i8)
         mid = lo + (hi-lo)/2_i8
         c = pcc(mid, mu, lambda)
         if (c >= p) then
            hi = mid
         else
            lo = mid
         end if
      end do
      q = hi
   end function qcc

   elemental function ieee_nan() result(x)
      use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function ieee_nan
end module ccd_distribution
