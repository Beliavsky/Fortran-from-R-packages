! SPDX-License-Identifier: GPL-2.0-or-later
module fer_special
   use fer_kinds, only : dp
   implicit none
   private
   public :: normal_pdf, normal_cdf, regularized_gamma_p, regularized_gamma_q
   public :: noncentral_chisq_cdf
contains
   elemental real(dp) function normal_pdf(x) result(y)
      real(dp), intent(in) :: x
      y = exp(-0.5_dp*x*x)/sqrt(2.0_dp*acos(-1.0_dp))
   end function normal_pdf

   elemental real(dp) function normal_cdf(x) result(y)
      real(dp), intent(in) :: x
      y = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   real(dp) function regularized_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      integer, parameter :: itmax = 10000
      real(dp), parameter :: eps = 4.0_dp*epsilon(1.0_dp)
      integer :: n
      real(dp) :: ap, del, sumv, b, c, d, h, an
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = 0.0_dp
      else if (x <= tiny(1.0_dp)) then
         p = 0.0_dp
      else if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do n = 1, itmax
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) <= abs(sumv)*eps) exit
         end do
         p = sumv*exp(-x + a*log(x) - log_gamma(a))
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
            if (abs(del-1.0_dp) <= eps) exit
         end do
         p = 1.0_dp - exp(-x + a*log(x) - log_gamma(a))*h
      end if
      p = max(0.0_dp, min(1.0_dp, p))
   end function regularized_gamma_p

   real(dp) function regularized_gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      q = 1.0_dp - regularized_gamma_p(a, x)
   end function regularized_gamma_q

   real(dp) function noncentral_chisq_cdf(x, df, lambda, lower_tail) result(cdf)
      real(dp), intent(in) :: x, df, lambda
      logical, intent(in), optional :: lower_tail
      integer, parameter :: maxit = 100000
      real(dp), parameter :: tol = 1.0e-14_dp
      logical :: lower
      integer :: j, j0
      real(dp) :: half_lam, w0, w, sumv, term
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      if (x <= 0.0_dp) then
         cdf = merge(0.0_dp, 1.0_dp, lower)
         return
      end if
      if (lambda <= 1.0e-14_dp) then
         cdf = regularized_gamma_p(0.5_dp*df, 0.5_dp*x)
         if (.not. lower) cdf = 1.0_dp-cdf
         return
      end if
      half_lam = 0.5_dp*lambda
      j0 = int(half_lam)
      w0 = exp(-half_lam + real(j0,dp)*log(half_lam) - log_gamma(real(j0+1,dp)))
      sumv = w0*regularized_gamma_p(0.5_dp*df + real(j0,dp), 0.5_dp*x)
      w = w0
      do j = j0-1, 0, -1
         w = w*real(j+1,dp)/half_lam
         term = w*regularized_gamma_p(0.5_dp*df + real(j,dp), 0.5_dp*x)
         sumv = sumv + term
      end do
      w = w0
      do j = j0+1, maxit
         w = w*half_lam/real(j,dp)
         term = w*regularized_gamma_p(0.5_dp*df + real(j,dp), 0.5_dp*x)
         sumv = sumv + term
         if (abs(term) < tol .and. w < tol) exit
      end do
      cdf = max(0.0_dp, min(1.0_dp, sumv))
      if (.not. lower) cdf = 1.0_dp-cdf
   end function noncentral_chisq_cdf
end module fer_special
