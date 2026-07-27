! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on SharpeR, copyright 2012-2025 Steven E. Pav.
module sharper_math
   use sharper_kinds, only: dp, pi, sqrt_two, huge_dp
   implicit none
   private

   public :: normal_pdf, normal_cdf, normal_quantile
   public :: gamma_p, gamma_q, beta_inc
   public :: chisq_pdf, chisq_cdf, chisq_quantile
   public :: ncchisq_cdf, ncchisq_quantile
   public :: student_t_pdf, student_t_cdf, student_t_quantile
   public :: nct_pdf, nct_cdf, nct_quantile
   public :: f_pdf, f_cdf, f_quantile
   public :: ncf_pdf, ncf_cdf, ncf_quantile
   public :: random_seed_scalar, random_normal, random_gamma
   public :: random_chisq, random_student_t, random_nct
   public :: random_f, random_ncf
   public :: binomial_coefficient, clamp_probability

contains

   pure elemental function clamp_probability(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = min(1.0_dp, max(0.0_dp, x))
   end function clamp_probability

   pure elemental function normal_pdf(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pi)
   end function normal_pdf

   pure elemental function normal_cdf(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = 0.5_dp * erfc(-x / sqrt_two)
   end function normal_cdf

   pure elemental function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: x, q, r
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
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow

      if (p <= 0.0_dp) then
         x = -huge_dp
      else if (p >= 1.0_dp) then
         x = huge_dp
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if

      if (p > 0.0_dp .and. p < 1.0_dp) then
         x = x - (normal_cdf(x)-p) / max(normal_pdf(x), tiny(1.0_dp))
      end if
   end function normal_quantile

   pure function gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      real(dp) :: p
      integer, parameter :: itmax = 10000
      real(dp), parameter :: eps = 2.0e-15_dp
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
      integer :: n
      real(dp) :: ap, del, sumv, b, c, d, h, an

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = 0.0_dp
         return
      end if
      if (abs(x) <= tiny(1.0_dp)) then
         p = 0.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do n = 1, itmax
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) <= abs(sumv)*eps) exit
         end do
         p = sumv*exp(-x+a*log(x)-log_gamma(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/fpmin
         d = 1.0_dp/b
         h = d
         do n = 1, itmax
            an = -real(n,dp)*(real(n,dp)-a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= eps) exit
         end do
         p = 1.0_dp - exp(-x+a*log(x)-log_gamma(a))*h
      end if
      p = clamp_probability(p)
   end function gamma_p

   pure function gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      real(dp) :: q
      q = clamp_probability(1.0_dp-gamma_p(a,x))
   end function gamma_q

   pure function beta_cf(a, b, x) result(cf)
      real(dp), intent(in) :: a, b, x
      real(dp) :: cf
      integer, parameter :: maxit = 10000
      real(dp), parameter :: eps = 2.0e-15_dp
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
      integer :: m, m2
      real(dp) :: aa, c, d, del, h, qab, qam, qap

      qab = a+b
      qap = a+1.0_dp
      qam = a-1.0_dp
      c = 1.0_dp
      d = 1.0_dp-qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d
      do m = 1, maxit
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x / &
              ((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c
         aa = -(a+real(m,dp))*(qab+real(m,dp))*x / &
              ((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) exit
      end do
      cf = h
   end function beta_cf

   pure function beta_inc(a, b, x) result(v)
      real(dp), intent(in) :: a, b, x
      real(dp) :: v, bt
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         v = 0.0_dp
      else if (x <= 0.0_dp) then
         v = 0.0_dp
      else if (x >= 1.0_dp) then
         v = 1.0_dp
      else
         bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b) + &
                  a*log(x)+b*log(1.0_dp-x))
         if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
            v = bt*beta_cf(a,b,x)/a
         else
            v = 1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
         end if
         v = clamp_probability(v)
      end if
   end function beta_inc

   pure elemental function chisq_pdf(x, df) result(v)
      real(dp), intent(in) :: x, df
      real(dp) :: v, a
      a = 0.5_dp*df
      if (x < 0.0_dp .or. df <= 0.0_dp) then
         v = 0.0_dp
      else if (abs(x) <= tiny(1.0_dp)) then
         if (a < 1.0_dp) then
            v = huge_dp
         else if (abs(a-1.0_dp) <= epsilon(1.0_dp)) then
            v = 0.5_dp
         else
            v = 0.0_dp
         end if
      else
         v = exp((a-1.0_dp)*log(x)-0.5_dp*x-a*log(2.0_dp)-log_gamma(a))
      end if
   end function chisq_pdf

   pure elemental function chisq_cdf(x, df) result(v)
      real(dp), intent(in) :: x, df
      real(dp) :: v
      if (x <= 0.0_dp) then
         v = 0.0_dp
      else
         v = gamma_p(0.5_dp*df,0.5_dp*x)
      end if
   end function chisq_cdf

   function chisq_quantile(p, df) result(x)
      real(dp), intent(in) :: p, df
      real(dp) :: x, lo, hi, mid
      integer :: iter
      if (p <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         x = huge_dp
         return
      end if
      lo = 0.0_dp
      hi = max(1.0_dp,df)
      do while (chisq_cdf(hi,df) < p .and. hi < huge_dp/4.0_dp)
         hi = 2.0_dp*hi
      end do
      do iter = 1, 160
         mid = 0.5_dp*(lo+hi)
         if (chisq_cdf(mid,df) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp*(lo+hi)
   end function chisq_quantile

   function ncchisq_cdf(x, df, lambda) result(v)
      real(dp), intent(in) :: x, df, lambda
      real(dp) :: v, lh, w, sumw
      integer :: j, jmax
      if (x <= 0.0_dp) then
         v = 0.0_dp
         return
      end if
      if (lambda <= 1.0e-14_dp) then
         v = chisq_cdf(x,df)
         return
      end if
      lh = 0.5_dp*lambda
      jmax = max(80,int(lh+14.0_dp*sqrt(lh+1.0_dp)+50.0_dp))
      v = 0.0_dp
      sumw = 0.0_dp
      do j = 0, jmax
         w = poisson_weight(j,lh)
         v = v+w*gamma_p(0.5_dp*df+real(j,dp),0.5_dp*x)
         sumw = sumw+w
      end do
      if (sumw > 0.0_dp) v = v/sumw
      v = clamp_probability(v)
   end function ncchisq_cdf

   function ncchisq_quantile(p, df, lambda) result(x)
      real(dp), intent(in) :: p, df, lambda
      real(dp) :: x, lo, hi, mid
      integer :: iter
      if (p <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         x = huge_dp
         return
      end if
      lo = 0.0_dp
      hi = max(1.0_dp,df+lambda)
      do while (ncchisq_cdf(hi,df,lambda) < p)
         hi = 2.0_dp*hi
      end do
      do iter = 1, 140
         mid = 0.5_dp*(lo+hi)
         if (ncchisq_cdf(mid,df,lambda) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp*(lo+hi)
   end function ncchisq_quantile

   pure elemental function student_t_pdf(x, df) result(v)
      real(dp), intent(in) :: x, df
      real(dp) :: v
      if (df <= 0.0_dp) then
         v = 0.0_dp
      else
         v = exp(log_gamma(0.5_dp*(df+1.0_dp))-log_gamma(0.5_dp*df) - &
                 0.5_dp*log(df*pi)-0.5_dp*(df+1.0_dp)*log(1.0_dp+x*x/df))
      end if
   end function student_t_pdf

   pure elemental function student_t_cdf(x, df) result(v)
      real(dp), intent(in) :: x, df
      real(dp) :: v, ib, xx
      if (df <= 0.0_dp) then
         v = 0.0_dp
      else if (abs(x) <= tiny(1.0_dp)) then
         v = 0.5_dp
      else
         xx = df/(df+x*x)
         ib = beta_inc(0.5_dp*df,0.5_dp,xx)
         if (x > 0.0_dp) then
            v = 1.0_dp-0.5_dp*ib
         else
            v = 0.5_dp*ib
         end if
      end if
      v = clamp_probability(v)
   end function student_t_cdf

   function student_t_quantile(p, df) result(x)
      real(dp), intent(in) :: p, df
      real(dp) :: x, lo, hi, mid
      integer :: iter
      if (p <= 0.0_dp) then
         x = -huge_dp
         return
      else if (p >= 1.0_dp) then
         x = huge_dp
         return
      end if
      lo = -1.0_dp
      hi = 1.0_dp
      do while (student_t_cdf(lo,df) > p)
         lo = 2.0_dp*lo
      end do
      do while (student_t_cdf(hi,df) < p)
         hi = 2.0_dp*hi
      end do
      do iter = 1, 140
         mid = 0.5_dp*(lo+hi)
         if (student_t_cdf(mid,df) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp*(lo+hi)
   end function student_t_quantile

   pure function nct_integrand(s, x, df, delta, density) result(v)
      real(dp), intent(in) :: s, x, df, delta
      logical, intent(in) :: density
      real(dp) :: v, y, a, z, logw
      y = exp(s)
      a = 0.5_dp*df
      z = x*sqrt(2.0_dp*y/df)-delta
      logw = a*s-y-log_gamma(a)
      if (density) then
         v = normal_pdf(z)*sqrt(2.0_dp*y/df)*exp(logw)
      else
         v = normal_cdf(z)*exp(logw)
      end if
   end function nct_integrand

   function nct_integral(x, df, delta, density) result(v)
      real(dp), intent(in) :: x, df, delta
      logical, intent(in) :: density
      real(dp) :: v, a, slo, shi, h, sumv, s
      integer, parameter :: nstep = 1200
      integer :: i
      if (df <= 0.0_dp) then
         v = 0.0_dp
         return
      end if
      a = 0.5_dp*df
      if (a > 4.0_dp) then
         slo = log(max(1.0e-12_dp,a-14.0_dp*sqrt(a)))
      else
         slo = -32.0_dp
      end if
      shi = log(a+14.0_dp*sqrt(a+1.0_dp)+80.0_dp)
      h = (shi-slo)/real(nstep,dp)
      sumv = nct_integrand(slo,x,df,delta,density) + &
             nct_integrand(shi,x,df,delta,density)
      do i = 1, nstep-1
         s = slo+real(i,dp)*h
         if (mod(i,2) == 0) then
            sumv = sumv+2.0_dp*nct_integrand(s,x,df,delta,density)
         else
            sumv = sumv+4.0_dp*nct_integrand(s,x,df,delta,density)
         end if
      end do
      v = sumv*h/3.0_dp
      if (.not. density) v = clamp_probability(v)
   end function nct_integral

   function nct_pdf(x, df, delta) result(v)
      real(dp), intent(in) :: x, df, delta
      real(dp) :: v
      if (abs(delta) < 1.0e-14_dp) then
         v = student_t_pdf(x,df)
      else
         v = max(0.0_dp,nct_integral(x,df,delta,.true.))
      end if
   end function nct_pdf

   function nct_cdf(x, df, delta) result(v)
      real(dp), intent(in) :: x, df, delta
      real(dp) :: v
      if (abs(delta) < 1.0e-14_dp) then
         v = student_t_cdf(x,df)
      else
         v = nct_integral(x,df,delta,.false.)
      end if
   end function nct_cdf

   function nct_quantile(p, df, delta) result(x)
      real(dp), intent(in) :: p, df, delta
      real(dp) :: x, lo, hi, mid
      integer :: iter
      if (p <= 0.0_dp) then
         x = -huge_dp
         return
      else if (p >= 1.0_dp) then
         x = huge_dp
         return
      end if
      lo = min(-1.0_dp,delta-8.0_dp)
      hi = max(1.0_dp,delta+8.0_dp)
      do while (nct_cdf(lo,df,delta) > p)
         lo = 2.0_dp*lo-1.0_dp
      end do
      do while (nct_cdf(hi,df,delta) < p)
         hi = 2.0_dp*hi+1.0_dp
      end do
      do iter = 1, 90
         mid = 0.5_dp*(lo+hi)
         if (nct_cdf(mid,df,delta) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp*(lo+hi)
   end function nct_quantile

   pure elemental function f_pdf(x, df1, df2) result(v)
      real(dp), intent(in) :: x, df1, df2
      real(dp) :: v, a, b
      if (x <= 0.0_dp .or. df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
         v = 0.0_dp
      else
         a = 0.5_dp*df1
         b = 0.5_dp*df2
         v = exp(a*log(df1/df2)+(a-1.0_dp)*log(x) - &
                 (a+b)*log(1.0_dp+df1*x/df2) - &
                 (log_gamma(a)+log_gamma(b)-log_gamma(a+b)))
      end if
   end function f_pdf

   pure elemental function f_cdf(x, df1, df2) result(v)
      real(dp), intent(in) :: x, df1, df2
      real(dp) :: v, z
      if (x <= 0.0_dp) then
         v = 0.0_dp
      else
         z = df1*x/(df1*x+df2)
         v = beta_inc(0.5_dp*df1,0.5_dp*df2,z)
      end if
   end function f_cdf

   function f_quantile(p, df1, df2) result(x)
      real(dp), intent(in) :: p, df1, df2
      real(dp) :: x, lo, hi, mid
      integer :: iter
      if (p <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         x = huge_dp
         return
      end if
      lo = 0.0_dp
      hi = 1.0_dp
      do while (f_cdf(hi,df1,df2) < p)
         hi = 2.0_dp*hi
      end do
      do iter = 1, 140
         mid = 0.5_dp*(lo+hi)
         if (f_cdf(mid,df1,df2) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp*(lo+hi)
   end function f_quantile

   pure function poisson_weight(j, lambda_half) result(w)
      integer, intent(in) :: j
      real(dp), intent(in) :: lambda_half
      real(dp) :: w
      if (lambda_half <= tiny(1.0_dp)) then
         if (j == 0) then
            w = 1.0_dp
         else
            w = 0.0_dp
         end if
      else
         w = exp(-lambda_half+real(j,dp)*log(lambda_half)-log_gamma(real(j+1,dp)))
      end if
   end function poisson_weight

   function ncf_cdf(x, df1, df2, lambda) result(v)
      real(dp), intent(in) :: x, df1, df2, lambda
      real(dp) :: v, z, lh, w, term, sumw
      integer :: j, jmax
      if (x <= 0.0_dp) then
         v = 0.0_dp
         return
      end if
      if (lambda <= 1.0e-14_dp) then
         v = f_cdf(x,df1,df2)
         return
      end if
      z = df1*x/(df1*x+df2)
      lh = 0.5_dp*lambda
      jmax = max(80,int(lh+14.0_dp*sqrt(lh+1.0_dp)+50.0_dp))
      v = 0.0_dp
      sumw = 0.0_dp
      do j = 0, jmax
         w = poisson_weight(j,lh)
         term = beta_inc(0.5_dp*df1+real(j,dp),0.5_dp*df2,z)
         v = v+w*term
         sumw = sumw+w
      end do
      if (sumw > 0.0_dp) v = v/sumw
      v = clamp_probability(v)
   end function ncf_cdf

   function ncf_pdf(x, df1, df2, lambda) result(v)
      real(dp), intent(in) :: x, df1, df2, lambda
      real(dp) :: v, z, dzdx, lh, w, a, b, logbeta, term, sumw
      integer :: j, jmax
      if (x <= 0.0_dp) then
         v = 0.0_dp
         return
      end if
      if (lambda <= 1.0e-14_dp) then
         v = f_pdf(x,df1,df2)
         return
      end if
      z = df1*x/(df1*x+df2)
      dzdx = df1*df2/(df1*x+df2)**2
      lh = 0.5_dp*lambda
      b = 0.5_dp*df2
      jmax = max(80,int(lh+14.0_dp*sqrt(lh+1.0_dp)+50.0_dp))
      v = 0.0_dp
      sumw = 0.0_dp
      do j = 0, jmax
         w = poisson_weight(j,lh)
         a = 0.5_dp*df1+real(j,dp)
         logbeta = log_gamma(a)+log_gamma(b)-log_gamma(a+b)
         term = exp((a-1.0_dp)*log(z)+(b-1.0_dp)*log(1.0_dp-z)-logbeta)*dzdx
         v = v+w*term
         sumw = sumw+w
      end do
      if (sumw > 0.0_dp) v = v/sumw
      v = max(0.0_dp,v)
   end function ncf_pdf

   function ncf_quantile(p, df1, df2, lambda) result(x)
      real(dp), intent(in) :: p, df1, df2, lambda
      real(dp) :: x, lo, hi, mid
      integer :: iter
      if (p <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         x = huge_dp
         return
      end if
      lo = 0.0_dp
      hi = max(1.0_dp,(df1+lambda)/df1)
      do while (ncf_cdf(hi,df1,df2,lambda) < p)
         hi = 2.0_dp*hi
      end do
      do iter = 1, 140
         mid = 0.5_dp*(lo+hi)
         if (ncf_cdf(mid,df1,df2,lambda) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp*(lo+hi)
   end function ncf_quantile

   subroutine random_seed_scalar(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: values(:)
      call random_seed(size=n)
      allocate(values(n))
      do i = 1, n
         values(i) = modulo(seed+104729*i,huge(1)-1)
         if (values(i) <= 0) values(i) = i
      end do
      call random_seed(put=values)
   end subroutine random_seed_scalar

   function random_normal() result(z)
      real(dp) :: z, u1, u2
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1,tiny(1.0_dp))
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function random_normal

   recursive function random_gamma(shape, scale) result(x)
      real(dp), intent(in) :: shape, scale
      real(dp) :: x, d, c, z, u, g
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         x = 0.0_dp
      else if (shape < 1.0_dp) then
         call random_number(u)
         g = random_gamma(shape+1.0_dp,1.0_dp)
         x = scale*g*u**(1.0_dp/shape)
      else
         d = shape-1.0_dp/3.0_dp
         c = 1.0_dp/sqrt(9.0_dp*d)
         do
            z = random_normal()
            if (1.0_dp+c*z <= 0.0_dp) cycle
            x = (1.0_dp+c*z)**3
            call random_number(u)
            if (u < 1.0_dp-0.0331_dp*z**4) exit
            if (log(u) < 0.5_dp*z*z+d*(1.0_dp-x+log(x))) exit
         end do
         x = scale*d*x
      end if
   end function random_gamma

   function random_chisq(df) result(x)
      real(dp), intent(in) :: df
      real(dp) :: x
      x = random_gamma(0.5_dp*df,2.0_dp)
   end function random_chisq

   function random_student_t(df) result(x)
      real(dp), intent(in) :: df
      real(dp) :: x
      x = random_normal()/sqrt(random_chisq(df)/df)
   end function random_student_t

   function random_nct(df, delta) result(x)
      real(dp), intent(in) :: df, delta
      real(dp) :: x
      x = (random_normal()+delta)/sqrt(random_chisq(df)/df)
   end function random_nct

   function random_f(df1, df2) result(x)
      real(dp), intent(in) :: df1, df2
      real(dp) :: x
      x = (random_chisq(df1)/df1)/(random_chisq(df2)/df2)
   end function random_f

   function random_ncf(df1, df2, lambda) result(x)
      real(dp), intent(in) :: df1, df2, lambda
      real(dp) :: x, z
      z = random_normal()+sqrt(max(0.0_dp,lambda))
      x = ((z*z+random_chisq(max(df1-1.0_dp,0.0_dp)))/df1) / &
          (random_chisq(df2)/df2)
   end function random_ncf

   pure function binomial_coefficient(n, k) result(v)
      integer, intent(in) :: n, k
      real(dp) :: v
      integer :: i, kk
      if (k < 0 .or. k > n) then
         v = 0.0_dp
         return
      end if
      kk = min(k,n-k)
      v = 1.0_dp
      do i = 1, kk
         v = v*real(n-kk+i,dp)/real(i,dp)
      end do
   end function binomial_coefficient

end module sharper_math
