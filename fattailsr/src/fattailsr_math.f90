! FatTailsR modern Fortran translation
! Copyright (C) 2014-2026 Patrice Kiener
! Licensed under GPL-2.0-only. See COPYING.
module fattailsr_math
   use fattailsr_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   private

   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: sqrt3 = sqrt(3.0_dp)
   real(dp), parameter, public :: kiener_scale = sqrt3/pi

   public :: logit, invlogit, kashp, ashp, dkashp_dx
   public :: beta_fn, incomplete_beta, regularized_beta
   public :: normal_pdf, normal_cdf, normal_quantile
   public :: clamp_probability, binomial_real

contains

   elemental pure function clamp_probability(p) result(q)
      real(dp), intent(in) :: p
      real(dp) :: q
      q = min(max(p, epsilon(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
   end function clamp_probability

   elemental pure function logit(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: x, q
      q = clamp_probability(p)
      x = log(q) - log(1.0_dp - q)
   end function logit

   elemental pure function invlogit(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: p
      if (x >= 0.0_dp) then
         p = 1.0_dp/(1.0_dp + exp(-x))
      else
         p = exp(x)/(1.0_dp + exp(x))
      end if
   end function invlogit

   elemental pure function ashp(x, k) result(z)
      real(dp), intent(in) :: x, k
      real(dp) :: z
      if (abs(k) <= tiny(1.0_dp)) then
         z = 0.0_dp
      else
         z = asinh(x/k)
      end if
   end function ashp

   elemental pure function kashp(x, k) result(z)
      real(dp), intent(in) :: x, k
      real(dp) :: z
      if (abs(k) <= tiny(1.0_dp)) then
         z = 0.0_dp
      else
         z = k*asinh(x/k)
      end if
   end function kashp

   elemental pure function dkashp_dx(x, k) result(z)
      real(dp), intent(in) :: x, k
      real(dp) :: z
      if (abs(k) <= tiny(1.0_dp)) then
         z = 0.0_dp
      else
         z = abs(k)/sqrt(x*x + k*k)
      end if
   end function dkashp_dx

   elemental pure function beta_fn(a, b) result(v)
      real(dp), intent(in) :: a, b
      real(dp) :: v
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         v = ieee_value(v, ieee_quiet_nan)
      else
         v = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
      end if
   end function beta_fn

   pure function beta_continued_fraction(a, b, x) result(cf)
      real(dp), intent(in) :: a, b, x
      real(dp) :: cf
      integer, parameter :: max_iter = 400
      real(dp), parameter :: tiny = 1.0e-300_dp
      real(dp), parameter :: tol = 4.0_dp*epsilon(1.0_dp)
      integer :: m, m2
      real(dp) :: aa, c, d, del, h, qab, qam, qap

      qab = a + b
      qap = a + 1.0_dp
      qam = a - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab*x/qap
      if (abs(d) < tiny) d = tiny
      d = 1.0_dp/d
      h = d

      do m = 1, max_iter
         m2 = 2*m
         aa = real(m, dp)*(b - real(m, dp))*x/&
              ((qam + real(m2, dp))*(a + real(m2, dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < tiny) d = tiny
         c = 1.0_dp + aa/c
         if (abs(c) < tiny) c = tiny
         d = 1.0_dp/d
         h = h*d*c

         aa = -(a + real(m, dp))*(qab + real(m, dp))*x/&
              ((a + real(m2, dp))*(qap + real(m2, dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < tiny) d = tiny
         c = 1.0_dp + aa/c
         if (abs(c) < tiny) c = tiny
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del - 1.0_dp) <= tol) exit
      end do
      cf = h
   end function beta_continued_fraction

   pure function regularized_beta(x, a, b) result(v)
      real(dp), intent(in) :: x, a, b
      real(dp) :: v, bt

      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         v = ieee_value(v, ieee_quiet_nan)
      else if (x <= 0.0_dp) then
         v = 0.0_dp
      else if (x >= 1.0_dp) then
         v = 1.0_dp
      else
         bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
                  a*log(x) + b*log(1.0_dp - x))
         if (x < (a + 1.0_dp)/(a + b + 2.0_dp)) then
            v = bt*beta_continued_fraction(a, b, x)/a
         else
            v = 1.0_dp - bt*beta_continued_fraction(b, a, 1.0_dp - x)/b
         end if
         v = min(max(v, 0.0_dp), 1.0_dp)
      end if
   end function regularized_beta

   pure function incomplete_beta(x, a, b) result(v)
      real(dp), intent(in) :: x, a, b
      real(dp) :: v
      v = regularized_beta(x, a, b)*beta_fn(a, b)
   end function incomplete_beta

   elemental pure function normal_pdf(x, mean, sd) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: mean, sd
      real(dp) :: v, mu, sigma, z
      mu = 0.0_dp
      sigma = 1.0_dp
      if (present(mean)) mu = mean
      if (present(sd)) sigma = sd
      z = (x - mu)/sigma
      v = exp(-0.5_dp*z*z)/(sqrt(2.0_dp*pi)*sigma)
   end function normal_pdf

   elemental pure function normal_cdf(x, mean, sd) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: mean, sd
      real(dp) :: v, mu, sigma
      mu = 0.0_dp
      sigma = 1.0_dp
      if (present(mean)) mu = mean
      if (present(sd)) sigma = sd
      v = 0.5_dp*erfc(-(x - mu)/(sigma*sqrt(2.0_dp)))
   end function normal_cdf

   pure function normal_quantile(p, mean, sd) result(x)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: mean, sd
      real(dp) :: x, q, r, mu, sigma
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
          4.374664141464968_dp, 2.938163982698783_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
          2.445134137142996_dp, 3.754408661907416_dp ]
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow

      mu = 0.0_dp
      sigma = 1.0_dp
      if (present(mean)) mu = mean
      if (present(sd)) sigma = sd
      q = clamp_probability(p)
      if (q < plow) then
         r = sqrt(-2.0_dp*log(q))
         x = (((((c(1)*r + c(2))*r + c(3))*r + c(4))*r + c(5))*r + c(6))/&
             ((((d(1)*r + d(2))*r + d(3))*r + d(4))*r + 1.0_dp)
      else if (q <= phigh) then
         r = q - 0.5_dp
         q = r*r
         x = (((((a(1)*q + a(2))*q + a(3))*q + a(4))*q + a(5))*q + a(6))*r/&
             (((((b(1)*q + b(2))*q + b(3))*q + b(4))*q + b(5))*q + 1.0_dp)
      else
         r = sqrt(-2.0_dp*log(1.0_dp - q))
         x = -(((((c(1)*r + c(2))*r + c(3))*r + c(4))*r + c(5))*r + c(6))/&
              ((((d(1)*r + d(2))*r + d(3))*r + d(4))*r + 1.0_dp)
      end if
      x = mu + sigma*x
   end function normal_quantile

   pure function binomial_real(n, k) result(v)
      integer, intent(in) :: n, k
      real(dp) :: v
      integer :: i, kk
      if (k < 0 .or. k > n) then
         v = 0.0_dp
         return
      end if
      kk = min(k, n-k)
      v = 1.0_dp
      do i = 1, kk
         v = v*real(n - kk + i, dp)/real(i, dp)
      end do
   end function binomial_real

end module fattailsr_math
