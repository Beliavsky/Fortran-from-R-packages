! SPDX-License-Identifier: GPL-2.0-only
module ppcor_special
   use ppcor_kinds, only : dp
   implicit none
   private
   public :: normal_cdf, student_t_cdf

contains

   pure real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   real(dp) function student_t_cdf(x, df) result(p)
      real(dp), intent(in) :: x, df
      real(dp) :: z, ib

      if (df <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      if (.not. (x > 0.0_dp .or. x < 0.0_dp)) then
         p = 0.5_dp
         return
      end if

      z = df / (df + x*x)
      ib = regularized_beta(z, 0.5_dp*df, 0.5_dp)
      if (x > 0.0_dp) then
         p = 1.0_dp - 0.5_dp*ib
      else
         p = 0.5_dp*ib
      end if
      p = min(1.0_dp, max(0.0_dp, p))
   end function student_t_cdf

   real(dp) function regularized_beta(x, a, b) result(value)
      real(dp), intent(in) :: x, a, b
      real(dp) :: front

      if (x <= 0.0_dp) then
         value = 0.0_dp
         return
      else if (x >= 1.0_dp) then
         value = 1.0_dp
         return
      end if

      front = exp(log_gamma(a+b) - log_gamma(a) - log_gamma(b) + &
                  a*log(x) + b*log(1.0_dp-x))
      if (x < (a + 1.0_dp)/(a + b + 2.0_dp)) then
         value = front * beta_continued_fraction(x, a, b) / a
      else
         value = 1.0_dp - front * beta_continued_fraction(1.0_dp-x, b, a) / b
      end if
      value = min(1.0_dp, max(0.0_dp, value))
   end function regularized_beta

   real(dp) function beta_continued_fraction(x, a, b) result(cf)
      real(dp), intent(in) :: x, a, b
      integer, parameter :: max_iter = 400
      real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
      integer :: m, m2
      real(dp) :: aa, c, d, delta, h, qab, qam, qap

      qab = a + b
      qap = a + 1.0_dp
      qam = a - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d

      do m = 1, max_iter
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x / &
              ((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c

         aa = -(a+real(m,dp))*(qab+real(m,dp))*x / &
              ((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         delta = d*c
         h = h*delta
         if (abs(delta-1.0_dp) <= eps) exit
      end do
      cf = h
   end function beta_continued_fraction

end module ppcor_special
