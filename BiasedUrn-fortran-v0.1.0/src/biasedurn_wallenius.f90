! BiasedUrn-fortran
! Wallenius' noncentral hypergeometric distribution.
! Upstream BiasedUrn copyright (c) 2002-2024 Agner Fog.
! License: GNU General Public License version 3.
module biasedurn_wallenius
   use biasedurn_kinds, only : dp
   use biasedurn_math, only : log_choose, wallenius_log_integral, quiet_nan, &
      positive_inf, rand_uniform
   use biasedurn_fisher, only : minhypergeo, maxhypergeo
   implicit none
   private

   public :: dwnchypergeo, pwnchypergeo, qwnchypergeo, rwnchypergeo
   public :: meanwnchypergeo, varwnchypergeo, modewnchypergeo
   public :: oddswnchypergeo, numwnchypergeo

contains

   function dwnchypergeo(x, m1, m2, n, odds, precision) result(p)
      integer, intent(in) :: x, m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(in), optional :: precision
      real(dp) :: p, logi, logp
      integer :: xmin, xmax
      integer :: xv(2), mv(2)
      real(dp) :: wv(2)

      p = 0.0_dp
      if (m1 < 0 .or. m2 < 0 .or. n < 0 .or. n > m1 + m2 .or. odds < 0.0_dp) then
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
         if (n > m2) then
            p = quiet_nan()
         else
            p = merge(1.0_dp, 0.0_dp, x == 0)
         end if
         return
      end if
      if (odds == 1.0_dp) then
         p = exp(log_choose(m1, x) + log_choose(m2, n - x) &
            - log_choose(m1 + m2, n))
         return
      end if
      if (n == m1 + m2) then
         p = merge(1.0_dp, 0.0_dp, x == m1)
         return
      end if

      xv = [x, n - x]
      mv = [m1, m2]
      wv = [odds, 1.0_dp]
      logi = wallenius_log_integral(xv, mv, wv, precision)
      if (logi /= logi) then
         p = quiet_nan()
         return
      end if
      logp = log_choose(m1, x) + log_choose(m2, n - x) + logi
      if (logp < log(tiny(1.0_dp))) then
         p = 0.0_dp
      else
         p = exp(logp)
      end if
      p = min(1.0_dp, max(0.0_dp, p))
   end function dwnchypergeo

   function pwnchypergeo(x, m1, m2, n, odds, precision, lower_tail) result(cdf)
      integer, intent(in) :: x, m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(in), optional :: precision
      logical, intent(in), optional :: lower_tail
      real(dp) :: cdf
      logical :: lower
      integer :: k, xmin, xmax

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
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
            cdf = cdf + dwnchypergeo(k, m1, m2, n, odds, precision)
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
            cdf = cdf + dwnchypergeo(k, m1, m2, n, odds, precision)
         end do
      end if
      cdf = min(1.0_dp, max(0.0_dp, cdf))
   end function pwnchypergeo

   integer function qwnchypergeo(prob, m1, m2, n, odds, precision, lower_tail) result(q)
      real(dp), intent(in) :: prob, odds
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in), optional :: precision
      logical, intent(in), optional :: lower_tail
      logical :: lower
      real(dp) :: c
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
         c = 0.0_dp
         do k = xmin, xmax
            c = c + dwnchypergeo(k, m1, m2, n, odds, precision)
            if (c >= prob) then
               q = k
               return
            end if
         end do
      else
         c = 1.0_dp
         do k = xmin, xmax
            c = c - dwnchypergeo(k, m1, m2, n, odds, precision)
            if (c <= prob) then
               q = k
               return
            end if
         end do
      end if
      q = xmax
   end function qwnchypergeo

   function rwnchypergeo(nran, m1, m2, n, odds, precision) result(draws)
      ! Exact urn-process simulation defining Wallenius' distribution.
      integer, intent(in) :: nran, m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(in), optional :: precision
      integer, allocatable :: draws(:)
      integer :: i, j, red, white, xr
      real(dp) :: wr, ww, u

      allocate(draws(max(0, nran)))
      if (nran <= 0) return
      do i = 1, nran
         red = m1
         white = m2
         xr = 0
         do j = 1, n
            wr = odds * real(red, dp)
            ww = real(white, dp)
            if (wr + ww <= 0.0_dp) then
               xr = -huge(1)
               exit
            end if
            u = rand_uniform() * (wr + ww)
            if (u < wr) then
               xr = xr + 1
               red = red - 1
            else
               white = white - 1
            end if
         end do
         draws(i) = xr
      end do
   end function rwnchypergeo

   subroutine wallenius_moments(m1, m2, n, odds, mean, variance, precision)
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
         p = dwnchypergeo(k, m1, m2, n, odds, precision)
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
   end subroutine wallenius_moments

   function meanwnchypergeo(m1, m2, n, odds, precision) result(mu)
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(in), optional :: precision
      real(dp) :: mu, v
      call wallenius_moments(m1, m2, n, odds, mu, v, precision)
   end function meanwnchypergeo

   function varwnchypergeo(m1, m2, n, odds, precision) result(v)
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in) :: odds
      real(dp), intent(in), optional :: precision
      real(dp) :: mu, v
      call wallenius_moments(m1, m2, n, odds, mu, v, precision)
   end function varwnchypergeo

   integer function modewnchypergeo(m1, m2, n, odds, precision) result(mode)
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
         p = dwnchypergeo(k, m1, m2, n, odds, precision)
         if (p > pbest) then
            pbest = p
            mode = k
         end if
      end do
   end function modewnchypergeo

   function oddswnchypergeo(mu, m1, m2, n, precision) result(odds)
      ! Manly approximation used by BiasedUrn.
      real(dp), intent(in) :: mu
      integer, intent(in) :: m1, m2, n
      real(dp), intent(in), optional :: precision
      real(dp) :: odds, a, b
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
         a = 1.0_dp - mu / real(m1, dp)
         b = 1.0_dp - (real(n, dp) - mu) / real(m2, dp)
         odds = log(a) / log(b)
      end if
   end function oddswnchypergeo

   function numwnchypergeo(mu, n, n_total, odds, precision) result(m)
      ! Solve Manly's mean approximation for urn composition.
      real(dp), intent(in) :: mu, odds
      integer, intent(in) :: n, n_total
      real(dp), intent(in), optional :: precision
      real(dp) :: m(2)
      real(dp) :: m1, m2, lastm1, mu2, z, zd
      integer :: iter

      if (n == 0 .or. odds < 0.0_dp) then
         m = quiet_nan()
         return
      end if
      if (n == n_total) then
         m = [mu, real(n_total, dp) - mu]
         return
      end if
      if (mu < 0.0_dp .or. mu > real(n, dp)) then
         m = quiet_nan()
         return
      end if
      if (mu == 0.0_dp) then
         m = [0.0_dp, real(n_total, dp)]
         return
      else if (mu == real(n, dp)) then
         m = [real(n_total, dp), 0.0_dp]
         return
      else if (odds == 0.0_dp) then
         m = quiet_nan()
         return
      end if

      mu2 = real(n, dp) - mu
      m1 = real(n_total, dp) * mu / real(n, dp)
      m2 = real(n_total, dp) - m1
      do iter = 1, 200
         lastm1 = m1
         z = log(1.0_dp - mu / m1) - odds * log(1.0_dp - mu2 / m2)
         zd = mu / (m1 * (m1 - mu)) &
            + odds * mu2 / (m2 * (m2 - mu2))
         m1 = m1 - z / zd
         if (m1 <= mu) m1 = 0.5_dp * (lastm1 + mu)
         m2 = real(n_total, dp) - m1
         if (m2 <= mu2) then
            m2 = 0.5_dp * (real(n_total, dp) - lastm1 + mu2)
            m1 = real(n_total, dp) - m2
         end if
         if (abs(m1 - lastm1) <= real(n_total, dp) * 1.0e-10_dp) exit
      end do
      m = [m1, real(n_total, dp) - m1]
   end function numwnchypergeo

end module biasedurn_wallenius
