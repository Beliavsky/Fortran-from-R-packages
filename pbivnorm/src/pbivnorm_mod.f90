! SPDX-License-Identifier: GPL-2.0-or-later
module pbivnorm_mod
   use, intrinsic :: iso_fortran_env, only : real64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   implicit none
   private

   integer, parameter, public :: dp = real64
   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp), parameter :: twopi = 2.0_dp*pi

   public :: pbivnorm
   public :: pbivnorm_recycle
   public :: bvn_upper_tail

contains

   elemental function pbivnorm(x, y, rho) result(prob)
      !! Standard bivariate normal CDF P(X <= x, Y <= y).
      real(dp), intent(in) :: x, y, rho
      real(dp) :: prob

      if (ieee_is_nan(x) .or. ieee_is_nan(y) .or. ieee_is_nan(rho)) then
         prob = x + y + rho
         return
      end if

      if (abs(rho) > 1.0_dp) then
         prob = ieee_value_like_nan()
         return
      end if

      if (.not. ieee_is_finite(x)) then
         if (x < 0.0_dp) then
            prob = 0.0_dp
         else
            prob = normal_cdf(y)
         end if
         return
      end if

      if (.not. ieee_is_finite(y)) then
         if (y < 0.0_dp) then
            prob = 0.0_dp
         else
            prob = normal_cdf(x)
         end if
         return
      end if

      if (rho >= 1.0_dp) then
         prob = normal_cdf(min(x, y))
      else if (rho <= -1.0_dp) then
         prob = max(0.0_dp, normal_cdf(x) - normal_cdf(-y))
      else
         prob = bvn_upper_tail(-x, -y, rho)
      end if

      prob = min(1.0_dp, max(0.0_dp, prob))
   end function pbivnorm

   subroutine pbivnorm_recycle(x, y, rho, prob, status)
      !! R-style vector recycling for x, y, and rho.
      !!
      !! prob must have size max(size(x), size(y), size(rho)).
      !! status = 0 on success, nonzero for invalid array sizes or rho.
      real(dp), intent(in) :: x(:), y(:), rho(:)
      real(dp), intent(out) :: prob(:)
      integer, intent(out), optional :: status
      integer :: i, n, istat

      istat = 0
      if (size(x) == 0 .or. size(y) == 0 .or. size(rho) == 0) then
         istat = 1
      else
         n = max(size(x), size(y), size(rho))
         if (size(prob) /= n) then
            istat = 2
         else
            do i = 1, n
               prob(i) = pbivnorm(x(1 + modulo(i - 1, size(x))), &
                                  y(1 + modulo(i - 1, size(y))), &
                                  rho(1 + modulo(i - 1, size(rho))))
               if (ieee_is_nan(prob(i))) istat = 3
            end do
         end if
      end if

      if (present(status)) status = istat
   end subroutine pbivnorm_recycle

   elemental function bvn_upper_tail(sh, sk, rho) result(prob)
      !! P(X > sh, Y > sk) for a standard bivariate normal pair.
      !! Modern free-form translation of the Genz/Ge MVBVU routine.
      real(dp), intent(in) :: sh, sk, rho
      real(dp) :: prob
      real(dp), parameter :: x(10, 3) = reshape([ &
         -0.9324695142031522_dp, -0.6612093864662647_dp, -0.2386191860831970_dp, &
          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         -0.9815606342467191_dp, -0.9041172563704750_dp, -0.7699026741943050_dp, &
         -0.5873179542866171_dp, -0.3678314989981802_dp, -0.1252334085114692_dp, &
          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         -0.9931285991850949_dp, -0.9639719272779138_dp, -0.9122344282513259_dp, &
         -0.8391169718222188_dp, -0.7463319064601508_dp, -0.6360536807265150_dp, &
         -0.5108670019508271_dp, -0.3737060887154196_dp, -0.2277858511416451_dp, &
         -0.07652652113349733_dp ], [10, 3])
      real(dp), parameter :: w(10, 3) = reshape([ &
          0.1713244923791705_dp, 0.3607615730481384_dp, 0.4679139345726904_dp, &
          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
          0.04717533638651177_dp, 0.1069393259953183_dp, 0.1600783285433464_dp, &
          0.2031674267230659_dp, 0.2334925365383547_dp, 0.2491470458134029_dp, &
          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
          0.01761400713915212_dp, 0.04060142980038694_dp, 0.06267204833410906_dp, &
          0.08327674157670475_dp, 0.1019301198172404_dp, 0.1181945319615184_dp, &
          0.1316886384491766_dp, 0.1420961093183821_dp, 0.1491729864726037_dp, &
          0.1527533871307259_dp ], [10, 3])
      real(dp) :: asv, a, b, c, d, rs, xs
      real(dp) :: sn, asr, h, k, bs, hs, hk
      integer :: i, lg, ng

      if (abs(rho) < 0.3_dp) then
         ng = 1
         lg = 3
      else if (abs(rho) < 0.75_dp) then
         ng = 2
         lg = 6
      else
         ng = 3
         lg = 10
      end if

      h = sh
      k = sk
      hk = h*k
      prob = 0.0_dp

      if (abs(rho) < 0.925_dp) then
         hs = 0.5_dp*(h*h + k*k)
         asr = asin(rho)
         do i = 1, lg
            sn = sin(0.5_dp*asr*(x(i, ng) + 1.0_dp))
            prob = prob + w(i, ng)*exp((sn*hk - hs)/(1.0_dp - sn*sn))
            sn = sin(0.5_dp*asr*(-x(i, ng) + 1.0_dp))
            prob = prob + w(i, ng)*exp((sn*hk - hs)/(1.0_dp - sn*sn))
         end do
         prob = prob*asr/(2.0_dp*twopi) + normal_cdf(-h)*normal_cdf(-k)
         return
      end if

      if (rho < 0.0_dp) then
         k = -k
         hk = -hk
      end if

      if (abs(rho) < 1.0_dp) then
         asv = (1.0_dp - rho)*(1.0_dp + rho)
         a = sqrt(asv)
         bs = (h - k)**2
         c = (4.0_dp - hk)/8.0_dp
         d = (12.0_dp - hk)/16.0_dp
         prob = a*exp(-0.5_dp*(bs/asv + hk))* &
            (1.0_dp - c*(bs - asv)*(1.0_dp - d*bs/5.0_dp)/3.0_dp + &
             c*d*asv*asv/5.0_dp)

         if (hk > -160.0_dp) then
            b = sqrt(bs)
            prob = prob - exp(-0.5_dp*hk)*sqrt(twopi)*normal_cdf(-b/a)*b* &
               (1.0_dp - c*bs*(1.0_dp - d*bs/5.0_dp)/3.0_dp)
         end if

         a = 0.5_dp*a
         do i = 1, lg
            xs = (a*(x(i, ng) + 1.0_dp))**2
            rs = sqrt(1.0_dp - xs)
            prob = prob + a*w(i, ng)*( &
               exp(-bs/(2.0_dp*xs) - hk/(1.0_dp + rs))/rs - &
               exp(-0.5_dp*(bs/xs + hk))*(1.0_dp + c*xs*(1.0_dp + d*xs)))

            xs = asv*(-x(i, ng) + 1.0_dp)**2/4.0_dp
            rs = sqrt(1.0_dp - xs)
            prob = prob + a*w(i, ng)*exp(-0.5_dp*(bs/xs + hk))* &
               (exp(-hk*(1.0_dp - rs)/(2.0_dp*(1.0_dp + rs)))/rs - &
                (1.0_dp + c*xs*(1.0_dp + d*xs)))
         end do
         prob = -prob/twopi
      end if

      if (rho > 0.0_dp) prob = prob + normal_cdf(-max(h, k))
      if (rho < 0.0_dp) prob = -prob + max(0.0_dp, normal_cdf(-h) - normal_cdf(-k))
   end function bvn_upper_tail

   elemental function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: p
      p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   elemental function ieee_value_like_nan() result(x)
      use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function ieee_value_like_nan

end module pbivnorm_mod
