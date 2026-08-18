! coda-fortran: computational translation of the R package coda.
! Original coda license: GPL (>= 2). This translation is GPL-2.0-or-later.
module coda_math
   use coda_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: mean_vec, variance_vec, covariance_matrix, sample_covariance
   public :: quantile_type7, sort_real, normal_quantile, f_quantile
   public :: regularized_beta, symmetric_max_eigenvalue, cholesky_lower
   public :: cramer_cdf, is_finite

contains

   pure real(dp) function mean_vec(x) result(m)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         m = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         m = sum(x) / real(size(x), dp)
      end if
   end function mean_vec

   pure real(dp) function variance_vec(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x) < 2) then
         v = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      m = mean_vec(x)
      v = sum((x - m)**2) / real(size(x) - 1, dp)
   end function variance_vec

   function covariance_matrix(x) result(cov)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable :: cov(:,:)
      real(dp), allocatable :: mu(:)
      integer :: n, p, i, j

      n = size(x, 1)
      p = size(x, 2)
      allocate(cov(p,p), mu(p))
      do j = 1, p
         mu(j) = mean_vec(x(:,j))
      end do
      if (n < 2) then
         cov = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      do j = 1, p
         do i = 1, p
            cov(i,j) = sum((x(:,i) - mu(i)) * (x(:,j) - mu(j))) / real(n - 1, dp)
         end do
      end do
   end function covariance_matrix

   pure real(dp) function sample_covariance(x, y) result(c)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: mx, my
      if (size(x) /= size(y) .or. size(x) < 2) then
         c = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      mx = mean_vec(x)
      my = mean_vec(y)
      c = sum((x - mx) * (y - my)) / real(size(x) - 1, dp)
   end function sample_covariance

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_real

   function quantile_type7(x, p) result(q)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: p
      real(dp) :: q, h, frac
      real(dp), allocatable :: y(:)
      integer :: n, lo

      n = size(x)
      if (n == 0 .or. p < 0.0_dp .or. p > 1.0_dp) then
         q = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      y = x
      call sort_real(y)
      if (n == 1 .or. p == 0.0_dp) then
         q = y(1)
         return
      end if
      if (p == 1.0_dp) then
         q = y(n)
         return
      end if
      h = 1.0_dp + real(n - 1, dp) * p
      lo = floor(h)
      frac = h - real(lo, dp)
      q = y(lo) + frac * (y(lo + 1) - y(lo))
   end function quantile_type7

   pure real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a1 = -3.969683028665376d1, a2 = 2.209460984245205d2
      real(dp), parameter :: a3 = -2.759285104469687d2, a4 = 1.383577518672690d2
      real(dp), parameter :: a5 = -3.066479806614716d1, a6 = 2.506628277459239d0
      real(dp), parameter :: b1 = -5.447609879822406d1, b2 = 1.615858368580409d2
      real(dp), parameter :: b3 = -1.556989798598866d2, b4 = 6.680131188771972d1
      real(dp), parameter :: b5 = -1.328068155288572d1
      real(dp), parameter :: c1 = -7.784894002430293d-3, c2 = -3.223964580411365d-1
      real(dp), parameter :: c3 = -2.400758277161838d0, c4 = -2.549732539343734d0
      real(dp), parameter :: c5 = 4.374664141464968d0, c6 = 2.938163982698783d0
      real(dp), parameter :: d1 = 7.784695709041462d-3, d2 = 3.224671290700398d-1
      real(dp), parameter :: d3 = 2.445134137142996d0, d4 = 3.754408661907416d0
      real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
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
         r = q * q
         x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6) * q / &
             (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp - p))
         x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
              ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_dp)
      end if
   end function normal_quantile

   pure real(dp) function beta_cont_frac(a, b, x) result(cf)
      real(dp), intent(in) :: a, b, x
      integer, parameter :: maxit = 300
      real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
      real(dp) :: qab, qap, qam, c, d, h, aa, del
      integer :: m, m2

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
         if (abs(del - 1.0_dp) < eps) exit
      end do
      cf = h
   end function beta_cont_frac

   pure real(dp) function regularized_beta(x, a, b) result(v)
      real(dp), intent(in) :: x, a, b
      real(dp) :: bt
      if (x <= 0.0_dp) then
         v = 0.0_dp
         return
      else if (x >= 1.0_dp) then
         v = 1.0_dp
         return
      end if
      bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * log(x) + b * log(1.0_dp - x))
      if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
         v = bt * beta_cont_frac(a, b, x) / a
      else
         v = 1.0_dp - bt * beta_cont_frac(b, a, 1.0_dp - x) / b
      end if
      v = max(0.0_dp, min(1.0_dp, v))
   end function regularized_beta

   pure real(dp) function f_cdf(f, d1, d2) result(p)
      real(dp), intent(in) :: f, d1, d2
      real(dp) :: x
      if (f <= 0.0_dp) then
         p = 0.0_dp
      else
         x = d1 * f / (d1 * f + d2)
         p = regularized_beta(x, 0.5_dp * d1, 0.5_dp * d2)
      end if
   end function f_cdf

   pure real(dp) function f_quantile(p, d1, d2) result(q)
      real(dp), intent(in) :: p, d1, d2
      real(dp) :: lo, hi, mid
      integer :: it
      if (p <= 0.0_dp) then
         q = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         q = huge(1.0_dp)
         return
      end if
      lo = 0.0_dp
      hi = 1.0_dp
      do while (f_cdf(hi, d1, d2) < p .and. hi < 1.0e100_dp)
         hi = hi * 2.0_dp
      end do
      do it = 1, 120
         mid = 0.5_dp * (lo + hi)
         if (f_cdf(mid, d1, d2) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      q = 0.5_dp * (lo + hi)
   end function f_quantile

   subroutine cholesky_lower(a, l, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: l(size(a,1),size(a,2))
      integer, intent(out) :: info
      integer :: n, i, j, k
      real(dp) :: s

      n = size(a,1)
      if (size(a,2) /= n) error stop "cholesky_lower: nonsquare matrix"
      l = 0.0_dp
      info = 0
      do i = 1, n
         do j = 1, i
            s = a(i,j)
            do k = 1, j - 1
               s = s - l(i,k) * l(j,k)
            end do
            if (i == j) then
               if (s <= 0.0_dp) then
                  info = i
                  return
               end if
               l(i,j) = sqrt(s)
            else
               l(i,j) = s / l(j,j)
            end if
         end do
      end do
   end subroutine cholesky_lower

   subroutine solve_lower(l, b, x)
      real(dp), intent(in) :: l(:,:), b(:,:)
      real(dp), intent(out) :: x(size(b,1),size(b,2))
      integer :: i, j, k, n
      real(dp) :: s
      n = size(l,1)
      x = 0.0_dp
      do j = 1, size(b,2)
         do i = 1, n
            s = b(i,j)
            do k = 1, i - 1
               s = s - l(i,k) * x(k,j)
            end do
            x(i,j) = s / l(i,i)
         end do
      end do
   end subroutine solve_lower

   subroutine solve_upper_from_lower(l, b, x)
      real(dp), intent(in) :: l(:,:), b(:,:)
      real(dp), intent(out) :: x(size(b,1),size(b,2))
      integer :: i, j, k, n
      real(dp) :: s
      n = size(l,1)
      x = 0.0_dp
      do j = 1, size(b,2)
         do i = n, 1, -1
            s = b(i,j)
            do k = i + 1, n
               s = s - l(k,i) * x(k,j)
            end do
            x(i,j) = s / l(i,i)
         end do
      end do
   end subroutine solve_upper_from_lower

   real(dp) function symmetric_max_eigenvalue(a) result(emax)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable :: m(:,:)
      integer :: n, p, q, i, j, iter
      real(dp) :: app, aqq, apq, phi, c, s, mip, miq, offmax

      n = size(a,1)
      m = 0.5_dp * (a + transpose(a))
      if (n == 1) then
         emax = m(1,1)
         return
      end if
      do iter = 1, 100 * n * n
         offmax = 0.0_dp
         p = 1
         q = 2
         do j = 2, n
            do i = 1, j - 1
               if (abs(m(i,j)) > offmax) then
                  offmax = abs(m(i,j))
                  p = i
                  q = j
               end if
            end do
         end do
         if (offmax < 1.0e-12_dp * max(1.0_dp, maxval(abs(m)))) exit
         app = m(p,p)
         aqq = m(q,q)
         apq = m(p,q)
         phi = 0.5_dp * atan2(2.0_dp * apq, aqq - app)
         c = cos(phi)
         s = sin(phi)
         do i = 1, n
            if (i /= p .and. i /= q) then
               mip = m(i,p)
               miq = m(i,q)
               m(i,p) = c * mip - s * miq
               m(p,i) = m(i,p)
               m(i,q) = s * mip + c * miq
               m(q,i) = m(i,q)
            end if
         end do
         m(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
         m(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
         m(p,q) = 0.0_dp
         m(q,p) = 0.0_dp
      end do
      emax = maxval([(m(i,i), i=1,n)])
   end function symmetric_max_eigenvalue

   real(dp) function modified_bessel_k_quarter(x) result(kv)
      real(dp), intent(in) :: x
      real(dp) :: ivp, ivm, term, sumv, xx, mu, corr
      integer :: k
      if (x <= 0.0_dp) then
         kv = huge(1.0_dp)
         return
      end if
      if (x > 2.0_dp) then
         mu = 0.25_dp
         corr = 1.0_dp + (4.0_dp*mu*mu - 1.0_dp)/(8.0_dp*x) + &
                (4.0_dp*mu*mu - 1.0_dp)*(4.0_dp*mu*mu - 9.0_dp)/(2.0_dp*(8.0_dp*x)**2) + &
                (4.0_dp*mu*mu - 1.0_dp)*(4.0_dp*mu*mu - 9.0_dp)*(4.0_dp*mu*mu - 25.0_dp)/ &
                (6.0_dp*(8.0_dp*x)**3)
         kv = sqrt(pi/(2.0_dp*x)) * exp(-x) * corr
         return
      end if
      xx = 0.25_dp * x * x
      term = exp(0.25_dp * log(0.5_dp*x) - log_gamma(1.25_dp))
      sumv = term
      do k = 1, 500
         term = term * xx / (real(k,dp) * (real(k,dp) + 0.25_dp))
         sumv = sumv + term
         if (abs(term) < 1.0e-16_dp * abs(sumv)) exit
      end do
      ivp = sumv
      term = exp(-0.25_dp * log(0.5_dp*x) - log_gamma(0.75_dp))
      sumv = term
      do k = 1, 500
         term = term * xx / (real(k,dp) * (real(k,dp) - 0.25_dp))
         sumv = sumv + term
         if (abs(term) < 1.0e-16_dp * abs(sumv)) exit
      end do
      ivm = sumv
      kv = pi * (ivm - ivp) / (2.0_dp * sin(0.25_dp*pi))
   end function modified_bessel_k_quarter

   real(dp) function cramer_cdf(q, eps) result(p)
      real(dp), intent(in) :: q
      real(dp), intent(in), optional :: eps
      real(dp) :: tol, logeps, z, u, term
      integer :: k
      if (q <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      tol = 1.0e-5_dp
      if (present(eps)) tol = eps
      logeps = log(tol)
      p = 0.0_dp
      do k = 0, 3
         z = gamma(real(k,dp) + 0.5_dp) * sqrt(real(4*k + 1,dp)) / &
             (gamma(real(k + 1,dp)) * pi**1.5_dp * sqrt(q))
         u = real((4*k + 1)**2, dp) / (16.0_dp * q)
         if (u <= -logeps) then
            term = z * exp(-u) * modified_bessel_k_quarter(u)
            p = p + term
         end if
      end do
      p = max(0.0_dp, min(1.0_dp, p))
   end function cramer_cdf

   pure logical function is_finite(x) result(ok)
      real(dp), intent(in) :: x
      ok = ieee_is_finite(x)
   end function is_finite

end module coda_math
