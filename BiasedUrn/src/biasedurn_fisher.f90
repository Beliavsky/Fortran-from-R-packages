! BiasedUrn-fortran
! Fisher's noncentral hypergeometric distribution.
! Upstream BiasedUrn copyright (c) 2002-2024 Agner Fog.
! License: GNU General Public License version 3.
module biasedurn_fisher
   use biasedurn_kinds, only : dp
   use biasedurn_math, only : log_choose, log_add, log_zero, quiet_nan, &
      positive_inf, rand_uniform
   implicit none
   private

   public :: dfnchypergeo, pfnchypergeo, qfnchypergeo, rfnchypergeo
   public :: meanfnchypergeo, varfnchypergeo, modefnchypergeo
   public :: oddsfnchypergeo, numfnchypergeo
   public :: minhypergeo, maxhypergeo

contains

   pure integer function minhypergeo(m1, m2, n) result(xmin)
      integer, intent(in) :: m1, m2, n
      xmin = max(n - m2, 0)
   end function minhypergeo

   pure integer function maxhypergeo(m1, m2, n) result(xmax)
      integer, intent(in) :: m1, m2, n
      xmax = min(m1, n)
   end function maxhypergeo

   function fisher_log_normalizer(m1, m2, n, odds) result(logz)
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in) :: odds
      real(dp) :: logz, t
      integer :: x, xmin, xmax

      xmin = minhypergeo(m1, m2, n)
      xmax = maxhypergeo(m1, m2, n)
      logz = log_zero
      do x = xmin, xmax
         if (odds == 0.0_dp .and. x > 0) cycle
         t = log_choose(m1, x) + log_choose(m2, n - x)
         if (x > 0) t = t + real(x, dp) * log(odds)
         logz = log_add(logz, t)
      end do
   end function fisher_log_normalizer

   function dfnchypergeo(x, m1, m2, n, odds, precision) result(p)
      integer, intent(in) :: x, m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(in), optional :: precision
      real(dp) :: p, lognum, logz
      integer :: xmin, xmax

      p = 0.0_dp
      if (.not. valid_parameters(m1, m2, n, odds)) then
         p = quiet_nan()
         return
      end if
      xmin = minhypergeo(m1, m2, n)
      xmax = maxhypergeo(m1, m2, n)
      if (x < xmin .or. x > xmax) return
      if (xmin == xmax) then
         p = 1.0_dp
         return
      end if
      if (odds == 0.0_dp) then
         p = merge(1.0_dp, 0.0_dp, x == xmin .and. xmin == 0)
         return
      end if
      if (odds == 1.0_dp) then
         p = exp(log_choose(m1, x) + log_choose(m2, n - x) &
            - log_choose(m1 + m2, n))
         return
      end if
      logz = fisher_log_normalizer(m1, m2, n, odds)
      lognum = log_choose(m1, x) + log_choose(m2, n - x) &
         + real(x, dp) * log(odds)
      p = exp(lognum - logz)
   end function dfnchypergeo

   function pfnchypergeo(x, m1, m2, n, odds, precision, lower_tail) result(cdf)
      integer, intent(in) :: x, m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(in), optional :: precision
      logical, intent(in), optional :: lower_tail
      real(dp) :: cdf
      logical :: lower
      integer :: k, xmin, xmax

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      if (.not. valid_parameters(m1, m2, n, odds)) then
         cdf = quiet_nan()
         return
      end if
      xmin = minhypergeo(m1, m2, n)
      xmax = maxhypergeo(m1, m2, n)
      if (lower) then
         if (x < xmin) then
            cdf = 0.0_dp
            return
         else if (x >= xmax) then
            cdf = 1.0_dp
            return
         end if
         cdf = 0.0_dp
         do k = xmin, x
            cdf = cdf + dfnchypergeo(k, m1, m2, n, odds, precision)
         end do
      else
         if (x < xmin) then
            cdf = 1.0_dp
            return
         else if (x >= xmax) then
            cdf = 0.0_dp
            return
         end if
         cdf = 0.0_dp
         do k = x + 1, xmax
            cdf = cdf + dfnchypergeo(k, m1, m2, n, odds, precision)
         end do
      end if
      cdf = min(1.0_dp, max(0.0_dp, cdf))
   end function pfnchypergeo

   integer function qfnchypergeo(prob, m1, m2, n, odds, precision, lower_tail) result(q)
      real(dp), intent(in) :: prob, odds
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in), optional :: precision
      logical, intent(in), optional :: lower_tail
      logical :: lower
      real(dp) :: target, c
      integer :: k, xmin, xmax

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      xmin = minhypergeo(m1, m2, n)
      xmax = maxhypergeo(m1, m2, n)
      if (prob <= 0.0_dp) then
         q = merge(xmin, xmax, lower)
         return
      else if (prob >= 1.0_dp) then
         q = merge(xmax, xmin, lower)
         return
      end if
      if (lower) then
         target = prob
         c = 0.0_dp
         do k = xmin, xmax
            c = c + dfnchypergeo(k, m1, m2, n, odds, precision)
            if (c >= target) then
               q = k
               return
            end if
         end do
         q = xmax
      else
         ! Lowest x such that P(X>x) <= prob.
         c = 1.0_dp
         do k = xmin, xmax
            c = c - dfnchypergeo(k, m1, m2, n, odds, precision)
            if (c <= prob) then
               q = k
               return
            end if
         end do
         q = xmax
      end if
   end function qfnchypergeo

   function rfnchypergeo(nran, m1, m2, n, odds, precision) result(draws)
      integer, intent(in) :: nran, m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(in), optional :: precision
      integer, allocatable :: draws(:)
      integer :: i, k, xmin, xmax
      real(dp) :: u, c

      allocate(draws(max(0, nran)))
      if (nran <= 0) return
      xmin = minhypergeo(m1, m2, n)
      xmax = maxhypergeo(m1, m2, n)
      do i = 1, nran
         u = rand_uniform()
         c = 0.0_dp
         draws(i) = xmax
         do k = xmin, xmax
            c = c + dfnchypergeo(k, m1, m2, n, odds, precision)
            if (u <= c) then
               draws(i) = k
               exit
            end if
         end do
      end do
   end function rfnchypergeo

   subroutine fisher_moments(m1, m2, n, odds, mean, variance, precision)
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(out) :: mean, variance
      real(dp), intent(in), optional :: precision
      integer :: k, xmin, xmax
      real(dp) :: p, s, sx, sx2
      xmin = minhypergeo(m1, m2, n)
      xmax = maxhypergeo(m1, m2, n)
      s = 0.0_dp
      sx = 0.0_dp
      sx2 = 0.0_dp
      do k = xmin, xmax
         p = dfnchypergeo(k, m1, m2, n, odds, precision)
         s = s + p
         sx = sx + real(k, dp) * p
         sx2 = sx2 + real(k * k, dp) * p
      end do
      if (s <= 0.0_dp) then
         mean = quiet_nan()
         variance = quiet_nan()
      else
         mean = sx / s
         variance = max(0.0_dp, sx2 / s - mean * mean)
      end if
   end subroutine fisher_moments

   function meanfnchypergeo(m1, m2, n, odds, precision) result(mu)
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(in), optional :: precision
      real(dp) :: mu, v
      call fisher_moments(m1, m2, n, odds, mu, v, precision)
   end function meanfnchypergeo

   function varfnchypergeo(m1, m2, n, odds, precision) result(v)
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(in), optional :: precision
      real(dp) :: mu, v
      call fisher_moments(m1, m2, n, odds, mu, v, precision)
   end function varfnchypergeo

   integer function modefnchypergeo(m1, m2, n, odds, precision) result(mode)
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(in), optional :: precision
      integer :: k, xmin, xmax
      real(dp) :: p, pbest
      xmin = minhypergeo(m1, m2, n)
      xmax = maxhypergeo(m1, m2, n)
      mode = xmin
      pbest = -1.0_dp
      do k = xmin, xmax
         p = dfnchypergeo(k, m1, m2, n, odds, precision)
         if (p > pbest) then
            pbest = p
            mode = k
         end if
      end do
   end function modefnchypergeo

   function oddsfnchypergeo(mu, m1, m2, n, precision) result(odds)
      real(dp), intent(in) :: mu
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in), optional :: precision
      real(dp) :: odds
      integer :: xmin, xmax
      xmin = minhypergeo(m1, m2, n)
      xmax = maxhypergeo(m1, m2, n)
      if (xmin == xmax) then
         odds = quiet_nan()
      else if (mu < real(xmin, dp) .or. mu > real(xmax, dp)) then
         odds = quiet_nan()
      else if (mu == real(xmin, dp)) then
         odds = 0.0_dp
      else if (mu == real(xmax, dp)) then
         odds = positive_inf()
      else
         odds = mu * (real(m2 - n, dp) + mu) / &
            ((real(m1, dp) - mu) * (real(n, dp) - mu))
      end if
   end function oddsfnchypergeo

   function numfnchypergeo(mu, n, n_total, odds, precision) result(m)
      real(dp), intent(in) :: mu, odds
      integer, intent(in) :: n, n_total
      real(dp), intent(in), optional :: precision
      real(dp) :: m(2), mu2, mu_o
      if (n == 0 .or. odds < 0.0_dp) then
         m = quiet_nan()
      else if (n == n_total) then
         m = [mu, real(n_total, dp) - mu]
      else if (mu < 0.0_dp .or. mu > real(n, dp)) then
         m = quiet_nan()
      else if (mu == 0.0_dp) then
         m = [0.0_dp, real(n_total, dp)]
      else if (mu == real(n, dp)) then
         m = [real(n_total, dp), 0.0_dp]
      else if (odds == 0.0_dp) then
         m = quiet_nan()
      else
         mu2 = real(n, dp) - mu
         mu_o = mu / odds
         m(1) = (mu_o * (real(n_total, dp) - mu2) + mu * mu2) / &
            (mu_o + mu2)
         m(2) = real(n_total, dp) - m(1)
      end if
   end function numfnchypergeo

   pure logical function valid_parameters(m1, m2, n, odds) result(ok)
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in) :: odds
      ok = m1 >= 0 .and. m2 >= 0 .and. n >= 0 .and. n <= m1 + m2 &
         .and. odds >= 0.0_dp .and. (odds > 0.0_dp .or. n <= m2)
   end function valid_parameters

end module biasedurn_fisher
