module lmtest_distributions
   use lmtest_kinds, only : dp
   implicit none
   private
   public :: normal_cdf, normal_sf, normal_quantile
   public :: chi_square_cdf, chi_square_sf
   public :: f_cdf, f_sf
   public :: student_t_cdf, student_t_sf, student_t_quantile

contains

   pure real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   pure real(dp) function normal_sf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp * erfc(x / sqrt(2.0_dp))
   end function normal_sf

   pure real(dp) function log_beta(a, b) result(v)
      real(dp), intent(in) :: a, b
      v = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
   end function log_beta

   real(dp) function beta_cont_frac(a, b, x) result(cf)
      real(dp), intent(in) :: a, b, x
      integer, parameter :: maxit = 400
      real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
      integer :: m, m2
      real(dp) :: aa, c, d, del, h, qab, qam, qap

      qab = a + b
      qap = a + 1.0_dp
      qam = a - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab * x / qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp / d
      h = d
      do m = 1, maxit
         m2 = 2 * m
         aa = real(m, dp) * (b - real(m, dp)) * x / &
              ((qam + real(m2, dp)) * (a + real(m2, dp)))
         d = 1.0_dp + aa * d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa / c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp / d
         h = h * d * c

         aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
              ((a + real(m2, dp)) * (qap + real(m2, dp)))
         d = 1.0_dp + aa * d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa / c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp / d
         del = d * c
         h = h * del
         if (abs(del - 1.0_dp) <= eps) exit
      end do
      cf = h
   end function beta_cont_frac

   real(dp) function regularized_beta(x, a, b) result(p)
      real(dp), intent(in) :: x, a, b
      real(dp) :: bt
      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      else if (x >= 1.0_dp) then
         p = 1.0_dp
         return
      end if
      bt = exp(a * log(x) + b * log(1.0_dp-x) - log_beta(a, b))
      if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
         p = bt * beta_cont_frac(a, b, x) / a
      else
         p = 1.0_dp - bt * beta_cont_frac(b, a, 1.0_dp - x) / b
      end if
      p = max(0.0_dp, min(1.0_dp, p))
   end function regularized_beta

   real(dp) function regularized_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      integer, parameter :: maxit = 10000
      real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
      integer :: n
      real(dp) :: ap, del, sumv, b, c, d, h, an

      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp / a
         del = sumv
         do n = 1, maxit
            ap = ap + 1.0_dp
            del = del * x / ap
            sumv = sumv + del
            if (abs(del) < abs(sumv) * eps) exit
         end do
         p = sumv * exp(-x + a * log(x) - log_gamma(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp / fpmin
         d = 1.0_dp / b
         h = d
         do n = 1, maxit
            an = -real(n, dp) * (real(n, dp) - a)
            b = b + 2.0_dp
            d = an * d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an / c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp / d
            del = d * c
            h = h * del
            if (abs(del - 1.0_dp) <= eps) exit
         end do
         p = 1.0_dp - exp(-x + a * log(x) - log_gamma(a)) * h
      end if
      p = max(0.0_dp, min(1.0_dp, p))
   end function regularized_gamma_p

   real(dp) function regularized_gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      integer, parameter :: maxit = 10000
      real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
      integer :: n
      real(dp) :: b, c, d, h, an, del

      if (x <= 0.0_dp) then
         q = 1.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         q = 1.0_dp - regularized_gamma_p(a, x)
         return
      end if
      b = x + 1.0_dp - a
      c = 1.0_dp / fpmin
      d = 1.0_dp / b
      h = d
      do n = 1, maxit
         an = -real(n, dp) * (real(n, dp) - a)
         b = b + 2.0_dp
         d = an * d + b
         if (abs(d) < fpmin) d = fpmin
         c = b + an / c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp / d
         del = d * c
         h = h * del
         if (abs(del - 1.0_dp) <= eps) exit
      end do
      q = exp(-x + a * log(x) - log_gamma(a)) * h
      q = max(0.0_dp, min(1.0_dp, q))
   end function regularized_gamma_q

   real(dp) function chi_square_cdf(x, df) result(p)
      real(dp), intent(in) :: x, df
      if (df <= 0.0_dp) then
         p = 0.0_dp
      else
         p = regularized_gamma_p(0.5_dp * df, 0.5_dp * max(x, 0.0_dp))
      end if
   end function chi_square_cdf

   real(dp) function chi_square_sf(x, df) result(p)
      real(dp), intent(in) :: x, df
      if (df <= 0.0_dp) then
         p = 1.0_dp
      else
         p = regularized_gamma_q(0.5_dp*df,0.5_dp*max(x,0.0_dp))
      end if
   end function chi_square_sf

   real(dp) function f_cdf(x, df1, df2) result(p)
      real(dp), intent(in) :: x, df1, df2
      real(dp) :: z
      if (x <= 0.0_dp .or. df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      z = df1 * x / (df1 * x + df2)
      p = regularized_beta(z, 0.5_dp * df1, 0.5_dp * df2)
   end function f_cdf

   real(dp) function f_sf(x, df1, df2) result(p)
      real(dp), intent(in) :: x, df1, df2
      real(dp) :: z
      if (x <= 0.0_dp .or. df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
         p = 1.0_dp
         return
      end if
      z = df2/(df2+df1*x)
      p = regularized_beta(z,0.5_dp*df2,0.5_dp*df1)
   end function f_sf

   real(dp) function student_t_cdf(t, df) result(p)
      real(dp), intent(in) :: t, df
      real(dp) :: x, ib
      if (df <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      if (abs(t) <= tiny(1.0_dp)) then
         p = 0.5_dp
         return
      end if
      x = df / (df + t * t)
      ib = regularized_beta(x, 0.5_dp * df, 0.5_dp)
      if (t > 0.0_dp) then
         p = 1.0_dp - 0.5_dp * ib
      else
         p = 0.5_dp * ib
      end if
   end function student_t_cdf

   real(dp) function student_t_sf(t, df) result(p)
      real(dp), intent(in) :: t, df
      real(dp) :: x, ib
      if (df <= 0.0_dp) then
         p = 1.0_dp
         return
      end if
      if (abs(t) <= tiny(1.0_dp)) then
         p = 0.5_dp
         return
      end if
      x = df/(df+t*t)
      ib = regularized_beta(x,0.5_dp*df,0.5_dp)
      if (t > 0.0_dp) then
         p = 0.5_dp*ib
      else
         p = 1.0_dp-0.5_dp*ib
      end if
   end function student_t_sf

   real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
          4.374664141464968_dp, 2.938163982698783_dp]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
          2.445134137142996_dp, 3.754408661907416_dp]
      real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
      real(dp) :: q, r

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp * log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q * q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
   end function normal_quantile

   real(dp) function student_t_quantile(p, df) result(x)
      real(dp), intent(in) :: p, df
      real(dp) :: lo, hi, mid
      integer :: i
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      end if
      lo = -1.0_dp
      hi = 1.0_dp
      do while (student_t_cdf(lo, df) > p)
         lo = 2.0_dp * lo
      end do
      do while (student_t_cdf(hi, df) < p)
         hi = 2.0_dp * hi
      end do
      do i = 1, 100
         mid = 0.5_dp * (lo + hi)
         if (student_t_cdf(mid, df) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp * (lo + hi)
   end function student_t_quantile

end module lmtest_distributions
