module psych_distributions
   use r_kinds, only : dp
   implicit none
   private
   public :: normal_cdf, normal_quantile, student_t_cdf, chi_square_cdf, f_cdf

contains

   pure elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x !! Standard-normal variate at which to evaluate the CDF.
      p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   pure elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p !! Probability in the open interval `(0,1)`.
      real(dp), parameter :: a1 = -3.969683028665376e1_dp
      real(dp), parameter :: a2 = 2.209460984245205e2_dp
      real(dp), parameter :: a3 = -2.759285104469687e2_dp
      real(dp), parameter :: a4 = 1.383577518672690e2_dp
      real(dp), parameter :: a5 = -3.066479806614716e1_dp
      real(dp), parameter :: a6 = 2.506628277459239_dp
      real(dp), parameter :: b1 = -5.447609879822406e1_dp
      real(dp), parameter :: b2 = 1.615858368580409e2_dp
      real(dp), parameter :: b3 = -1.556989798598866e2_dp
      real(dp), parameter :: b4 = 6.680131188771972e1_dp
      real(dp), parameter :: b5 = -1.328068155288572e1_dp
      real(dp), parameter :: c1 = -7.784894002430293e-3_dp
      real(dp), parameter :: c2 = -3.223964580411365e-1_dp
      real(dp), parameter :: c3 = -2.400758277161838_dp
      real(dp), parameter :: c4 = -2.549732539343734_dp
      real(dp), parameter :: c5 = 4.374664141464968_dp
      real(dp), parameter :: c6 = 2.938163982698783_dp
      real(dp), parameter :: d1 = 7.784695709041462e-3_dp
      real(dp), parameter :: d2 = 3.224671290700398e-1_dp
      real(dp), parameter :: d3 = 2.445134137142996_dp
      real(dp), parameter :: d4 = 3.754408661907416_dp
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow
      real(dp) :: q, r

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp * log(p))
         x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
            ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q*q
         x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q / &
            (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp - p))
         x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
            ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_dp)
      end if
   end function normal_quantile

   pure real(dp) function beta_continued_fraction(a, b, x) result(cf)
      real(dp), intent(in) :: a !! First beta-shape parameter, strictly positive.
      real(dp), intent(in) :: b !! Second beta-shape parameter, strictly positive.
      real(dp), intent(in) :: x !! Evaluation point in `[0,1]`.
      integer, parameter :: max_iter = 300
      real(dp), parameter :: eps = 3.0e-14_dp
      real(dp), parameter :: fpmin = 1.0e-300_dp
      real(dp) :: aa, c, d, del, h, qab, qam, qap
      integer :: m, m2

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
         aa = real(m, dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c
         aa = -(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) exit
      end do
      cf = h
   end function beta_continued_fraction

   pure real(dp) function regularized_beta(x, a, b) result(value)
      real(dp), intent(in) :: x !! Evaluation point in `[0,1]`.
      real(dp), intent(in) :: a !! First beta-shape parameter, strictly positive.
      real(dp), intent(in) :: b !! Second beta-shape parameter, strictly positive.
      real(dp) :: bt

      if (x <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (x >= 1.0_dp) then
         value = 1.0_dp
         return
      end if
      bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         value = bt*beta_continued_fraction(a,b,x)/a
      else
         value = 1.0_dp - bt*beta_continued_fraction(b,a,1.0_dp-x)/b
      end if
   end function regularized_beta

   pure elemental real(dp) function student_t_cdf(t, df) result(p)
      real(dp), intent(in) :: t !! Student-t variate.
      real(dp), intent(in) :: df !! Positive degrees of freedom.
      real(dp) :: x, ib

      if (df <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      x = df/(df+t*t)
      ib = regularized_beta(x, 0.5_dp*df, 0.5_dp)
      if (t >= 0.0_dp) then
         p = 1.0_dp - 0.5_dp*ib
      else
         p = 0.5_dp*ib
      end if
   end function student_t_cdf

   pure real(dp) function regularized_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a !! Positive gamma-shape parameter.
      real(dp), intent(in) :: x !! Nonnegative evaluation point.
      integer, parameter :: max_iter = 500
      real(dp), parameter :: eps = 3.0e-14_dp
      real(dp), parameter :: fpmin = 1.0e-300_dp
      real(dp) :: ap, del, sumv, b, c, d, h, an
      integer :: i

      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do i = 1, max_iter
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) < abs(sumv)*eps) exit
         end do
         p = sumv*exp(-x+a*log(x)-log_gamma(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/fpmin
         d = 1.0_dp/b
         h = d
         do i = 1, max_iter
            an = -real(i,dp)*(real(i,dp)-a)
            b = b + 2.0_dp
            d = an*d+b
            if (abs(d) < fpmin) d = fpmin
            c = b+an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) < eps) exit
         end do
         p = 1.0_dp - exp(-x+a*log(x)-log_gamma(a))*h
      end if
      p = min(1.0_dp, max(0.0_dp, p))
   end function regularized_gamma_p

   pure elemental real(dp) function chi_square_cdf(x, df) result(p)
      real(dp), intent(in) :: x !! Nonnegative chi-square statistic.
      real(dp), intent(in) :: df !! Positive degrees of freedom.
      if (x <= 0.0_dp) then
         p = 0.0_dp
      else
         p = regularized_gamma_p(0.5_dp*df, 0.5_dp*x)
      end if
   end function chi_square_cdf

   pure elemental real(dp) function f_cdf(x, df1, df2) result(p)
      real(dp), intent(in) :: x !! Nonnegative F statistic.
      real(dp), intent(in) :: df1 !! Positive numerator degrees of freedom.
      real(dp), intent(in) :: df2 !! Positive denominator degrees of freedom.
      real(dp) :: z
      if (x <= 0.0_dp) then
         p = 0.0_dp
      else
         z = df1*x/(df1*x+df2)
         p = regularized_beta(z, 0.5_dp*df1, 0.5_dp*df2)
      end if
   end function f_cdf

end module psych_distributions
