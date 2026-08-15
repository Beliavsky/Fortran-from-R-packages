! distr-fortran -- computational translation of the R package distr.
! Copyright (C) 2005-2025 distr authors.
! SPDX-License-Identifier: LGPL-3.0-only
module distr_special
   use distr_kinds, only : dp, pi, sqrt2, sqrt2pi, eps_dp, nan_dp
   implicit none
   private
   public :: normal_pdf_std, normal_cdf_std, normal_quantile_std
   public :: regularized_gamma_p, regularized_gamma_q, regularized_beta
   public :: beta_log_density, gamma_log_density
   public :: digamma_value, inverse_digamma
   public :: log_choose, is_integer_value, clamp01
   public :: poisson_weighted_beta, poisson_weighted_chisq
   public :: noncentral_t_cdf, noncentral_t_density

contains

   elemental real(dp) function clamp01(x) result(y)
      real(dp), intent(in) :: x
      y = min(1.0_dp, max(0.0_dp, x))
   end function clamp01

   elemental logical function is_integer_value(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x - anint(x)) <= 16.0_dp*eps_dp*max(1.0_dp, abs(x))
   end function is_integer_value

   elemental real(dp) function normal_pdf_std(x) result(y)
      real(dp), intent(in) :: x
      y = exp(-0.5_dp*x*x)/sqrt2pi
   end function normal_pdf_std

   elemental real(dp) function normal_cdf_std(x) result(y)
      real(dp), intent(in) :: x
      y = 0.5_dp*erfc(-x/sqrt2)
   end function normal_cdf_std

   elemental real(dp) function normal_quantile_std(p) result(x)
      ! Peter J. Acklam's rational approximation, followed by one Halley step.
      real(dp), intent(in) :: p
      real(dp), parameter :: a1=-3.969683028665376e1_dp, a2=2.209460984245205e2_dp
      real(dp), parameter :: a3=-2.759285104469687e2_dp, a4=1.383577518672690e2_dp
      real(dp), parameter :: a5=-3.066479806614716e1_dp, a6=2.506628277459239_dp
      real(dp), parameter :: b1=-5.447609879822406e1_dp, b2=1.615858368580409e2_dp
      real(dp), parameter :: b3=-1.556989798598866e2_dp, b4=6.680131188771972e1_dp
      real(dp), parameter :: b5=-1.328068155288572e1_dp
      real(dp), parameter :: c1=-7.784894002430293e-3_dp, c2=-3.223964580411365e-1_dp
      real(dp), parameter :: c3=-2.400758277161838_dp, c4=-2.549732539343734_dp
      real(dp), parameter :: c5=4.374664141464968_dp, c6=2.938163982698783_dp
      real(dp), parameter :: d1=7.784695709041462e-3_dp, d2=3.224671290700398e-1_dp
      real(dp), parameter :: d3=2.445134137142996_dp, d4=3.754408661907416_dp
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      real(dp) :: q, r, e, u
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      end if
      if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
             ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q*q
         x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / &
             (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
              ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      end if
      e = normal_cdf_std(x) - p
      u = e*sqrt2pi*exp(0.5_dp*x*x)
      x = x - u/(1.0_dp + 0.5_dp*x*u)
   end function normal_quantile_std

   real(dp) function regularized_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      integer, parameter :: itmax = 100000
      real(dp), parameter :: tol = 8.0_dp*eps_dp
      integer :: n
      real(dp) :: ap, del, sumv, gln, b, c, d, h, an
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = nan_dp()
         return
      end if
      if (x == 0.0_dp) then
         p = 0.0_dp
         return
      end if
      gln = log_gamma(a)
      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do n = 1, itmax
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) <= abs(sumv)*tol) exit
         end do
         p = sumv*exp(-x + a*log(x) - gln)
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/tiny(1.0_dp)
         d = 1.0_dp/b
         h = d
         do n = 1, itmax
            an = -real(n,dp)*(real(n,dp)-a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
            c = b + an/c
            if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del - 1.0_dp) <= tol) exit
         end do
         p = 1.0_dp - exp(-x + a*log(x) - gln)*h
      end if
      p = clamp01(p)
   end function regularized_gamma_p

   real(dp) function regularized_gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      integer, parameter :: itmax = 100000
      real(dp), parameter :: tol = 8.0_dp*eps_dp
      integer :: n
      real(dp) :: b, c, d, h, an, del, gln
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         q = nan_dp()
         return
      end if
      if (x == 0.0_dp) then
         q = 1.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         q = 1.0_dp - regularized_gamma_p(a,x)
         q = clamp01(q)
         return
      end if
      gln = log_gamma(a)
      b = x + 1.0_dp - a
      c = 1.0_dp/tiny(1.0_dp)
      d = 1.0_dp/b
      h = d
      do n = 1, itmax
         an = -real(n,dp)*(real(n,dp)-a)
         b = b + 2.0_dp
         d = an*d + b
         if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
         c = b + an/c
         if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del - 1.0_dp) <= tol) exit
      end do
      q = exp(-x + a*log(x) - gln)*h
      q = clamp01(q)
   end function regularized_gamma_q

   real(dp) function beta_cf(a, b, x) result(cf)
      real(dp), intent(in) :: a, b, x
      integer, parameter :: maxit = 100000
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps_dp, tol = 8.0_dp*eps_dp
      integer :: m, m2
      real(dp) :: aa, c, d, del, h, qab, qam, qap
      qab = a+b
      qap = a+1.0_dp
      qam = a-1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d
      do m = 1, maxit
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x / ((qam+real(m2,dp))*(a+real(m2,dp)))
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
         del = d*c
         h = h*del
         if (abs(del - 1.0_dp) <= tol) exit
      end do
      cf = h
   end function beta_cf

   real(dp) function regularized_beta(x, a, b) result(p)
      real(dp), intent(in) :: x, a, b
      real(dp) :: bt
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         p = nan_dp()
         return
      end if
      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      else if (x >= 1.0_dp) then
         p = 1.0_dp
         return
      end if
      bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b) + a*log(x) + b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         p = bt*beta_cf(a,b,x)/a
      else
         p = 1.0_dp - bt*beta_cf(b,a,1.0_dp-x)/b
      end if
      p = clamp01(p)
   end function regularized_beta

   elemental real(dp) function beta_log_density(x, a, b) result(v)
      real(dp), intent(in) :: x, a, b
      if (x < 0.0_dp .or. x > 1.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp) then
         v = -huge(1.0_dp)
      else if (x == 0.0_dp) then
         if (a < 1.0_dp) then
            v = huge(1.0_dp)
         else if (a == 1.0_dp) then
            v = log_gamma(a+b)-log_gamma(a)-log_gamma(b)
         else
            v = -huge(1.0_dp)
         end if
      else if (x == 1.0_dp) then
         if (b < 1.0_dp) then
            v = huge(1.0_dp)
         else if (b == 1.0_dp) then
            v = log_gamma(a+b)-log_gamma(a)-log_gamma(b)
         else
            v = -huge(1.0_dp)
         end if
      else
         v = (a-1.0_dp)*log(x)+(b-1.0_dp)*log(1.0_dp-x) + &
             log_gamma(a+b)-log_gamma(a)-log_gamma(b)
      end if
   end function beta_log_density

   elemental real(dp) function gamma_log_density(x, shape, scale) result(v)
      real(dp), intent(in) :: x, shape, scale
      if (shape <= 0.0_dp .or. scale <= 0.0_dp .or. x < 0.0_dp) then
         v = -huge(1.0_dp)
      else if (x == 0.0_dp) then
         if (shape < 1.0_dp) then
            v = huge(1.0_dp)
         else if (shape == 1.0_dp) then
            v = -log(scale)
         else
            v = -huge(1.0_dp)
         end if
      else
         v = (shape-1.0_dp)*log(x) - x/scale - log_gamma(shape) - shape*log(scale)
      end if
   end function gamma_log_density

   elemental real(dp) function log_choose(n, k) result(v)
      integer, intent(in) :: n, k
      if (k < 0 .or. k > n .or. n < 0) then
         v = -huge(1.0_dp)
      else
         v = log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp))
      end if
   end function log_choose

   recursive elemental real(dp) function digamma_value(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: z, r, inv, inv2
      if (x <= 0.0_dp) then
         if (abs(x-anint(x)) < 4.0_dp*eps_dp) then
            y = nan_dp()
            return
         end if
         y = digamma_value(1.0_dp-x) - pi/tan(pi*x)
         return
      end if
      z = x
      r = 0.0_dp
      do while (z < 8.0_dp)
         r = r - 1.0_dp/z
         z = z + 1.0_dp
      end do
      inv = 1.0_dp/z
      inv2 = inv*inv
      y = r + log(z) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - &
          inv2*(1.0_dp/120.0_dp - inv2*(1.0_dp/252.0_dp - &
          inv2*(1.0_dp/240.0_dp - inv2/132.0_dp))))
   end function digamma_value

   elemental real(dp) function trigamma_value(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: z, r, inv, inv2
      if (x <= 0.0_dp) then
         y = nan_dp()
         return
      end if
      z = x
      r = 0.0_dp
      do while (z < 8.0_dp)
         r = r + 1.0_dp/(z*z)
         z = z + 1.0_dp
      end do
      inv = 1.0_dp/z
      inv2 = inv*inv
      y = r + inv + 0.5_dp*inv2 + inv2*inv/6.0_dp - inv2*inv2*inv/30.0_dp + &
          inv2*inv2*inv2*inv/42.0_dp - inv2**4*inv/30.0_dp + &
          5.0_dp*inv2**5*inv/66.0_dp
   end function trigamma_value

   real(dp) function inverse_digamma(y) result(x)
      ! Inverse of digamma, replacing distr's interpolation table by Newton iteration.
      real(dp), intent(in) :: y
      integer :: it
      real(dp) :: step
      if (y >= -2.22_dp) then
         x = exp(y) + 0.5_dp
      else
         x = -1.0_dp/(y + 0.5772156649015328606_dp)
      end if
      x = max(x, tiny(1.0_dp)**0.25_dp)
      do it = 1, 30
         step = (digamma_value(x)-y)/trigamma_value(x)
         x = x - step
         if (x <= 0.0_dp) x = 0.5_dp*(x+step)
         if (abs(step) <= 16.0_dp*eps_dp*max(1.0_dp,x)) exit
      end do
   end function inverse_digamma

   real(dp) function poisson_mode_weight(lambda, k) result(w)
      real(dp), intent(in) :: lambda
      integer, intent(in) :: k
      if (lambda == 0.0_dp) then
         w = merge(1.0_dp, 0.0_dp, k == 0)
      else
         w = exp(-lambda + real(k,dp)*log(lambda) - log_gamma(real(k+1,dp)))
      end if
   end function poisson_mode_weight

   real(dp) function poisson_weighted_beta(x, a, b, ncp, density, upper_tail) result(v)
      real(dp), intent(in) :: x, a, b, ncp
      logical, intent(in) :: density
      logical, intent(in), optional :: upper_tail
      integer :: k0, k, nsmall
      real(dp) :: lam, w, term, sumv, wu, wd, contrib
      logical :: upper
      upper = .false.
      if (present(upper_tail)) upper = upper_tail
      if (ncp <= 0.0_dp) then
         if (density) then
            v = exp(beta_log_density(x,a,b))
         else
            if (upper) then
               v = regularized_beta(1.0_dp-x,b,a)
            else
               v = regularized_beta(x,a,b)
            end if
         end if
         return
      end if
      lam = 0.5_dp*ncp
      k0 = int(floor(lam))
      w = poisson_mode_weight(lam,k0)
      if (density) then
         term = exp(beta_log_density(x,a+real(k0,dp),b))
      else
         if (upper) then
            term = regularized_beta(1.0_dp-x,b,a+real(k0,dp))
         else
            term = regularized_beta(x,a+real(k0,dp),b)
         end if
      end if
      sumv = w*term
      wd = w
      nsmall = 0
      do k = k0-1, 0, -1
         wd = wd*real(k+1,dp)/lam
         if (density) then
            term = exp(beta_log_density(x,a+real(k,dp),b))
         else
            if (upper) then
               term = regularized_beta(1.0_dp-x,b,a+real(k,dp))
            else
               term = regularized_beta(x,a+real(k,dp),b)
            end if
         end if
         contrib = wd*term
         sumv = sumv + contrib
         if (wd < 1.0e-16_dp .and. abs(contrib) < 1.0e-15_dp*max(sumv,1.0e-300_dp)) exit
      end do
      wu = w
      do k = k0+1, k0+200000
         wu = wu*lam/real(k,dp)
         if (density) then
            term = exp(beta_log_density(x,a+real(k,dp),b))
         else
            if (upper) then
               term = regularized_beta(1.0_dp-x,b,a+real(k,dp))
            else
               term = regularized_beta(x,a+real(k,dp),b)
            end if
         end if
         contrib = wu*term
         sumv = sumv + contrib
         if (wu < 1.0e-16_dp .and. abs(contrib) < 1.0e-15_dp*max(sumv,1.0e-300_dp)) then
            nsmall = nsmall + 1
            if (nsmall >= 4) exit
         else
            nsmall = 0
         end if
         if (wu == 0.0_dp) exit
      end do
      v = sumv
      if (.not. density) v = clamp01(v)
   end function poisson_weighted_beta

   real(dp) function central_chisq_density(x, df) result(v)
      real(dp), intent(in) :: x, df
      if (x < 0.0_dp .or. df <= 0.0_dp) then
         v = 0.0_dp
      else
         v = exp(gamma_log_density(x,0.5_dp*df,2.0_dp))
      end if
   end function central_chisq_density

   real(dp) function poisson_weighted_chisq(x, df, ncp, density, upper_tail) result(v)
      real(dp), intent(in) :: x, df, ncp
      logical, intent(in) :: density
      logical, intent(in), optional :: upper_tail
      integer :: k0, k, nsmall
      real(dp) :: lam, w, term, sumv, wu, wd, contrib
      logical :: upper
      upper = .false.
      if (present(upper_tail)) upper = upper_tail
      if (x < 0.0_dp) then
         if (density) then
            v = 0.0_dp
         else
            v = merge(1.0_dp,0.0_dp,upper)
         end if
         return
      end if
      if (ncp <= 0.0_dp) then
         if (density) then
            v = central_chisq_density(x,df)
         else
            if (upper) then
               v = regularized_gamma_q(0.5_dp*df,0.5_dp*x)
            else
               v = regularized_gamma_p(0.5_dp*df,0.5_dp*x)
            end if
         end if
         return
      end if
      lam = 0.5_dp*ncp
      k0 = int(floor(lam))
      w = poisson_mode_weight(lam,k0)
      if (density) then
         term = central_chisq_density(x,df+2.0_dp*real(k0,dp))
      else
         if (upper) then
            term = regularized_gamma_q(0.5_dp*df+real(k0,dp),0.5_dp*x)
         else
            term = regularized_gamma_p(0.5_dp*df+real(k0,dp),0.5_dp*x)
         end if
      end if
      sumv = w*term
      wd = w
      do k = k0-1, 0, -1
         wd = wd*real(k+1,dp)/lam
         if (density) then
            term = central_chisq_density(x,df+2.0_dp*real(k,dp))
         else
            if (upper) then
               term = regularized_gamma_q(0.5_dp*df+real(k,dp),0.5_dp*x)
            else
               term = regularized_gamma_p(0.5_dp*df+real(k,dp),0.5_dp*x)
            end if
         end if
         contrib = wd*term
         sumv = sumv + contrib
         if (wd < 1.0e-16_dp .and. abs(contrib) < 1.0e-15_dp*max(sumv,1.0e-300_dp)) exit
      end do
      wu = w
      nsmall = 0
      do k = k0+1, k0+200000
         wu = wu*lam/real(k,dp)
         if (density) then
            term = central_chisq_density(x,df+2.0_dp*real(k,dp))
         else
            if (upper) then
               term = regularized_gamma_q(0.5_dp*df+real(k,dp),0.5_dp*x)
            else
               term = regularized_gamma_p(0.5_dp*df+real(k,dp),0.5_dp*x)
            end if
         end if
         contrib = wu*term
         sumv = sumv + contrib
         if (wu < 1.0e-16_dp .and. abs(contrib) < 1.0e-15_dp*max(sumv,1.0e-300_dp)) then
            nsmall = nsmall + 1
            if (nsmall >= 4) exit
         else
            nsmall = 0
         end if
         if (wu == 0.0_dp) exit
      end do
      v = sumv
      if (.not. density) v = clamp01(v)
   end function poisson_weighted_chisq

   real(dp) function noncentral_t_integral(t, df, ncp, density) result(v)
      real(dp), intent(in) :: t, df, ncp
      logical, intent(in) :: density
      integer, parameter :: ng = 96
      integer :: i, j
      real(dp) :: z, z1, p1, pp, u1, u2, w, vv, logf, arg, term, sumv
      real(dp), parameter :: tol = 2.0e-15_dp
      sumv = 0.0_dp
      do i = 1, (ng+1)/2
         z = cos(pi*(real(i,dp)-0.25_dp)/(real(ng,dp)+0.5_dp))
         ! Newton iteration for Legendre root.
         do j = 1, 20
            p1 = 1.0_dp
            pp = 0.0_dp
            call legendre_eval(ng,z,p1,pp)
            z1 = z
            z = z1 - p1/pp
            if (abs(z-z1) <= tol) exit
         end do
         call legendre_eval(ng,z,p1,pp)
         w = 2.0_dp/((1.0_dp-z*z)*pp*pp)
         u1 = 0.5_dp*(1.0_dp-z)
         u2 = 0.5_dp*(1.0_dp+z)
         sumv = sumv + 0.5_dp*w*(nct_u_integrand(u1,t,df,ncp,density) + &
                                  nct_u_integrand(u2,t,df,ncp,density))
      end do
      v = sumv
   contains
      subroutine legendre_eval(n,z,pn,dpn)
         integer, intent(in) :: n
         real(dp), intent(in) :: z
         real(dp), intent(out) :: pn, dpn
         integer :: k
         real(dp) :: p0, pcur, pnext
         p0 = 1.0_dp
         pcur = z
         if (n == 0) then
            pn = 1.0_dp; dpn = 0.0_dp; return
         else if (n == 1) then
            pn = z; dpn = 1.0_dp; return
         end if
         do k = 2, n
            pnext = ((2.0_dp*real(k,dp)-1.0_dp)*z*pcur-(real(k,dp)-1.0_dp)*p0)/real(k,dp)
            p0 = pcur
            pcur = pnext
         end do
         pn = pcur
         dpn = real(n,dp)*(z*pn-p0)/(z*z-1.0_dp)
      end subroutine legendre_eval

      real(dp) function nct_u_integrand(u,t0,nu,delta,want_d) result(g)
         real(dp), intent(in) :: u,t0,nu,delta
         logical, intent(in) :: want_d
         real(dp) :: xchi, lpdf, a0, jac
         if (u <= 0.0_dp .or. u >= 1.0_dp) then
            g = 0.0_dp; return
         end if
         xchi = u/(1.0_dp-u)
         jac = 1.0_dp/(1.0_dp-u)**2
         lpdf = gamma_log_density(xchi,0.5_dp*nu,2.0_dp)
         if (lpdf < log(tiny(1.0_dp))) then
            g = 0.0_dp; return
         end if
         a0 = t0*sqrt(xchi/nu)-delta
         if (want_d) then
            g = normal_pdf_std(a0)*sqrt(xchi/nu)*exp(lpdf)*jac
         else
            g = normal_cdf_std(a0)*exp(lpdf)*jac
         end if
      end function nct_u_integrand
   end function noncentral_t_integral

   real(dp) function noncentral_t_cdf(t, df, ncp) result(v)
      real(dp), intent(in) :: t, df, ncp
      if (df <= 0.0_dp) then
         v = nan_dp()
      else if (ncp == 0.0_dp) then
         if (t == 0.0_dp) then
            v = 0.5_dp
         else if (t > 0.0_dp) then
            v = 1.0_dp - 0.5_dp*regularized_beta(df/(df+t*t),0.5_dp*df,0.5_dp)
         else
            v = 0.5_dp*regularized_beta(df/(df+t*t),0.5_dp*df,0.5_dp)
         end if
      else
         v = clamp01(noncentral_t_integral(t,df,ncp,.false.))
      end if
   end function noncentral_t_cdf

   real(dp) function noncentral_t_density(t, df, ncp) result(v)
      real(dp), intent(in) :: t, df, ncp
      if (df <= 0.0_dp) then
         v = nan_dp()
      else if (ncp == 0.0_dp) then
         v = exp(log_gamma(0.5_dp*(df+1.0_dp))-log_gamma(0.5_dp*df) - &
             0.5_dp*log(df*pi) - 0.5_dp*(df+1.0_dp)*log(1.0_dp+t*t/df))
      else
         v = noncentral_t_integral(t,df,ncp,.true.)
      end if
   end function noncentral_t_density

end module distr_special
