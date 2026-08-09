! SPDX-License-Identifier: GPL-2.0-or-later
module coneproj_stats
   use coneproj_kinds, only : dp
   implicit none
   private
   public :: seed_rng, randn, regularized_beta, student_t_cdf, student_t_quantile

contains

   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729 * i, huge(1) - 1)
         if (put(i) == 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_rng

   real(dp) function randn() result(z)
      real(dp) :: u1, u2
      real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp * log(u1)) * cos(twopi * u2)
   end function randn

   real(dp) function regularized_beta(x, a, b) result(val)
      real(dp), intent(in) :: x, a, b
      real(dp) :: bt
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         val = 0.0_dp
         return
      end if
      if (x <= 0.0_dp) then
         val = 0.0_dp
         return
      end if
      if (x >= 1.0_dp) then
         val = 1.0_dp
         return
      end if
      bt = exp(log_gamma(a+b) - log_gamma(a) - log_gamma(b) + a*log(x) + b*log(1.0_dp-x))
      if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
         val = bt * beta_cf(a,b,x) / a
      else
         val = 1.0_dp - bt * beta_cf(b,a,1.0_dp-x) / b
      end if
      val = min(1.0_dp, max(0.0_dp, val))
   end function regularized_beta

   real(dp) function beta_cf(a, b, x) result(h)
      real(dp), intent(in) :: a, b, x
      integer, parameter :: maxit = 200
      real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
      real(dp) :: qab, qap, qam, c, d, aa, del
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
         m2 = 2*m
         aa = real(m,dp) * (b-real(m,dp)) * x / &
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
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) exit
      end do
   end function beta_cf

   real(dp) function student_t_cdf(t, df) result(p)
      real(dp), intent(in) :: t, df
      real(dp) :: x, ib
      if (df <= 0.0_dp) then
         p = 0.5_dp
         return
      end if
      if (abs(t) <= tiny(1.0_dp)) then
         p = 0.5_dp
         return
      end if
      x = df / (df + t*t)
      ib = regularized_beta(x, 0.5_dp*df, 0.5_dp)
      if (t > 0.0_dp) then
         p = 1.0_dp - 0.5_dp*ib
      else
         p = 0.5_dp*ib
      end if
   end function student_t_cdf

   real(dp) function student_t_quantile(p, df) result(q)
      real(dp), intent(in) :: p, df
      real(dp) :: lo, hi, mid
      integer :: i
      if (p <= 0.0_dp) then
         q = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         q = huge(1.0_dp)
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
      do i = 1, 120
         mid = 0.5_dp*(lo+hi)
         if (student_t_cdf(mid,df) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      q = 0.5_dp*(lo+hi)
   end function student_t_quantile

end module coneproj_stats
