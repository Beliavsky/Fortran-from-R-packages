! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
module poilog_distribution
   use poilog_kinds, only : dp
   use poilog_math, only : pi_dp, normal_pdf, poisson_pmf, safe_exp, is_finite_dp
   use poilog_quadrature, only : integrate_gk15
   implicit none
   private
   public :: dpoilog, dpoilog_vec, dbipoilog, dbipoilog_vec

contains

   function dpoilog(n, mu, sig) result(p)
      integer, intent(in) :: n
      real(dp), intent(in) :: mu, sig
      real(dp) :: p
      real(dp) :: p_log, p_norm

      if (n < 0 .or. .not. is_finite_dp(mu) .or. .not. is_finite_dp(sig) .or. sig <= 0.0_dp) then
         p = 0.0_dp
         return
      end if

      p_log = poilog_loglambda(n, mu, sig*sig)
      if (n < 8) then
         p_norm = poilog_normal(n, mu, sig)
         p = max(p_log, p_norm)
      else
         p = p_log
         if (.not. is_finite_dp(p) .or. p < 0.0_dp) p = poilog_normal(n, mu, sig)
      end if
      if (.not. is_finite_dp(p) .or. p < 0.0_dp) p = 0.0_dp
      p = min(p, 1.0_dp)
   end function dpoilog

   function dpoilog_vec(n, mu, sig) result(p)
      integer, intent(in) :: n(:)
      real(dp), intent(in) :: mu, sig
      real(dp) :: p(size(n))
      integer :: i
      do i = 1, size(n)
         p(i) = dpoilog(n(i), mu, sig)
      end do
   end function dpoilog_vec

   function dbipoilog(n1, n2, mu1, mu2, sig1, sig2, rho) result(p)
      integer, intent(in) :: n1, n2
      real(dp), intent(in) :: mu1, mu2, sig1, sig2, rho
      real(dp) :: p
      real(dp) :: var1, var2, cond_var, cond_sig, m, a, b, fac
      integer :: stat

      if (n1 < 0 .or. n2 < 0) then
         p = 0.0_dp
         return
      end if
      if (.not. all([is_finite_dp(mu1), is_finite_dp(mu2), is_finite_dp(sig1), &
                     is_finite_dp(sig2), is_finite_dp(rho)])) then
         p = 0.0_dp
         return
      end if
      if (sig1 <= 0.0_dp .or. sig2 <= 0.0_dp .or. abs(rho) > 1.0_dp) then
         p = 0.0_dp
         return
      end if

      var1 = sig1*sig1
      var2 = sig2*sig2
      if (1.0_dp-rho*rho <= 64.0_dp*epsilon(1.0_dp)) then
         p = bipoilog_singular(n1,n2,mu1,mu2,sig1,sig2,rho)
         return
      end if

      cond_var = var2*(1.0_dp-rho*rho)
      cond_sig = sqrt(cond_var)
      m = maxf(n1, mu1, var1)
      a = lower_bound(n1, m, mu1, var1)
      b = upper_bound(n1, m, mu1, var1)
      fac = log_gamma(real(n1+1,dp))

      p = integrate_gk15(f_outer, a, b, 1.0e-7_dp, 1.0e-7_dp, stat) / sqrt(2.0_dp*pi_dp*var1)
      if (stat /= 0 .or. .not. is_finite_dp(p) .or. p < 0.0_dp) then
         p = bipoilog_normal(n1,n2,mu1,mu2,sig1,sig2,rho)
      end if
      p = max(0.0_dp, min(p, 1.0_dp))

   contains
      function f_outer(z) result(v)
         real(dp), intent(in) :: z
         real(dp) :: v, cmu, logv, ez
         cmu = mu2 + rho*sqrt(var2/var1)*(z-mu1)
         ez = exp_for_penalty(z)
         if (ez >= huge(1.0_dp)/4.0_dp) then
            v = 0.0_dp
            return
         end if
         logv = real(n1,dp)*z - ez - fac - 0.5_dp*(z-mu1)*(z-mu1)/var1
         if (logv < log(tiny(1.0_dp))) then
            v = 0.0_dp
         else
            v = dpoilog(n2, cmu, cond_sig)*exp(logv)
         end if
      end function f_outer
   end function dbipoilog

   function dbipoilog_vec(n1, n2, mu1, mu2, sig1, sig2, rho) result(p)
      integer, intent(in) :: n1(:), n2(:)
      real(dp), intent(in) :: mu1, mu2, sig1, sig2, rho
      real(dp) :: p(size(n1))
      integer :: i
      if (size(n1) /= size(n2)) then
         p = 0.0_dp
         return
      end if
      do i = 1, size(n1)
         p(i) = dbipoilog(n1(i),n2(i),mu1,mu2,sig1,sig2,rho)
      end do
   end function dbipoilog_vec

   function poilog_normal(n, mu, sig) result(p)
      integer, intent(in) :: n
      real(dp), intent(in) :: mu, sig
      real(dp) :: p
      integer :: stat
      p = integrate_gk15(f, -12.0_dp, 12.0_dp, 1.0e-10_dp, 1.0e-9_dp, stat)
      if (stat /= 0 .or. .not. is_finite_dp(p)) p = max(0.0_dp,p)
   contains
      function f(u) result(v)
         real(dp), intent(in) :: u
         real(dp) :: v, eta, lambda
         eta = mu + sig*u
         if (eta >= log(huge(1.0_dp))) then
            v = 0.0_dp
         else
            lambda = safe_exp(eta)
            v = normal_pdf(u)*poisson_pmf(n,lambda)
         end if
      end function f
   end function poilog_normal

   function bipoilog_normal(n1,n2,mu1,mu2,sig1,sig2,rho) result(p)
      integer, intent(in) :: n1,n2
      real(dp), intent(in) :: mu1,mu2,sig1,sig2,rho
      real(dp) :: p
      real(dp) :: cond_sig
      integer :: stat
      cond_sig = sig2*sqrt(max(0.0_dp,1.0_dp-rho*rho))
      p = integrate_gk15(f, -12.0_dp, 12.0_dp, 1.0e-9_dp, 1.0e-8_dp, stat)
      if (.not. is_finite_dp(p) .or. p < 0.0_dp) p = 0.0_dp
   contains
      function f(u) result(v)
         real(dp), intent(in) :: u
         real(dp) :: v, cmu
         cmu = mu2 + sig2*rho*u
         v = normal_pdf(u)*poisson_pmf(n1,safe_exp(mu1+sig1*u))*dpoilog(n2,cmu,cond_sig)
      end function f
   end function bipoilog_normal

   function bipoilog_singular(n1,n2,mu1,mu2,sig1,sig2,rho) result(p)
      integer, intent(in) :: n1,n2
      real(dp), intent(in) :: mu1,mu2,sig1,sig2,rho
      real(dp) :: p
      integer :: stat
      p = integrate_gk15(f, -12.0_dp, 12.0_dp, 1.0e-10_dp, 1.0e-8_dp, stat)
      if (.not. is_finite_dp(p) .or. p < 0.0_dp) p = 0.0_dp
   contains
      function f(u) result(v)
         real(dp), intent(in) :: u
         real(dp) :: v
         v = normal_pdf(u)*poisson_pmf(n1,safe_exp(mu1+sig1*u))* &
             poisson_pmf(n2,safe_exp(mu2+sig2*sign(1.0_dp,rho)*u))
      end function f
   end function bipoilog_singular

   function poilog_loglambda(n, mu, var) result(p)
      integer, intent(in) :: n
      real(dp), intent(in) :: mu, var
      real(dp) :: p, m, a, b, fac
      integer :: stat
      m = maxf(n,mu,var)
      a = lower_bound(n,m,mu,var)
      b = upper_bound(n,m,mu,var)
      fac = log_gamma(real(n+1,dp))
      p = integrate_gk15(f, a, b, 1.0e-9_dp, 1.0e-8_dp, stat)/sqrt(2.0_dp*pi_dp*var)
      if (stat /= 0 .or. .not. is_finite_dp(p)) p = 0.0_dp
   contains
      function f(z) result(v)
         real(dp), intent(in) :: z
         real(dp) :: v, e, lv
         e = exp_for_penalty(z)
         if (e >= huge(1.0_dp)/4.0_dp) then
            v = 0.0_dp
            return
         end if
         lv = real(n,dp)*z - e - 0.5_dp*(z-mu)*(z-mu)/var - fac
         if (lv <= log(tiny(1.0_dp))) then
            v = 0.0_dp
         else
            v = exp(lv)
         end if
      end function f
   end function poilog_loglambda

   pure real(dp) function maxf(n, mu, var) result(z)
      integer, intent(in) :: n
      real(dp), intent(in) :: mu, var
      real(dp) :: d
      z = 0.0_dp
      d = 100.0_dp
      do while (d > 1.0e-5_dp)
         if (real(n-1,dp)-exp_for_penalty(z)-(z-mu)/var > 0.0_dp) then
            z = z+d
         else
            z = z-d
         end if
         d = 0.5_dp*d
      end do
   end function maxf

   pure real(dp) function upper_bound(n,m,mu,var) result(z)
      integer, intent(in) :: n
      real(dp), intent(in) :: m,mu,var
      real(dp) :: d,mf,g
      mf = real(n-1,dp)*m-exp_for_penalty(m)-0.5_dp*(m-mu)*(m-mu)/var
      z = m+20.0_dp
      d = 10.0_dp
      do while (d > 1.0e-6_dp)
         g = real(n-1,dp)*z-exp_for_penalty(z)-0.5_dp*(z-mu)*(z-mu)/var-mf+log(1.0e6_dp)
         if (g > 0.0_dp) then
            z = z+d
         else
            z = z-d
         end if
         d = 0.5_dp*d
      end do
   end function upper_bound

   pure real(dp) function lower_bound(n,m,mu,var) result(z)
      integer, intent(in) :: n
      real(dp), intent(in) :: m,mu,var
      real(dp) :: d,mf,g
      mf = real(n-1,dp)*m-exp_for_penalty(m)-0.5_dp*(m-mu)*(m-mu)/var
      z = m-20.0_dp
      d = 10.0_dp
      do while (d > 1.0e-6_dp)
         g = real(n-1,dp)*z-exp_for_penalty(z)-0.5_dp*(z-mu)*(z-mu)/var-mf+log(1.0e6_dp)
         if (g > 0.0_dp) then
            z = z-d
         else
            z = z+d
         end if
         d = 0.5_dp*d
      end do
   end function lower_bound

   pure real(dp) function exp_for_penalty(z) result(v)
      real(dp), intent(in) :: z
      if (z > log(huge(1.0_dp))-2.0_dp) then
         v = huge(1.0_dp)/2.0_dp
      else if (z < log(tiny(1.0_dp))) then
         v = 0.0_dp
      else
         v = exp(z)
      end if
   end function exp_for_penalty

end module poilog_distribution
