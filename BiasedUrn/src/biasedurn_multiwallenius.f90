! BiasedUrn-fortran
! Multivariate Wallenius' noncentral hypergeometric distribution.
! Upstream BiasedUrn copyright (c) 2002-2024 Agner Fog.
! License: GNU General Public License version 3.
module biasedurn_multiwallenius
   use biasedurn_kinds, only : dp
   use biasedurn_math, only : log_choose, wallenius_log_integral, quiet_nan, &
      positive_inf, rand_uniform
   implicit none
   private

   public :: dmwnchypergeo, rmwnchypergeo, momentsmwnchypergeo
   public :: meanmwnchypergeo, varmwnchypergeo
   public :: oddsmwnchypergeo, nummwnchypergeo

contains

   function dmwnchypergeo(x, m, n, odds, precision) result(p)
      integer, intent(in) :: x(:), m(:), n
      real(dp), intent(in) :: odds(:)
      real(dp), intent(in), optional :: precision
      real(dp) :: p, logi, logp
      integer :: i

      p = 0.0_dp
      if (.not. valid_multi(x, m, n, odds)) then
         if (size(x) /= size(m) .or. size(m) /= size(odds)) p = quiet_nan()
         return
      end if
      if (n == 0) then
         p = 1.0_dp
         return
      end if
      if (n == sum(m)) then
         if (all(x == m)) p = 1.0_dp
         return
      end if
      if (all(abs(odds - odds(1)) <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(odds(1))))) then
         logp = -log_choose(sum(m), n)
         do i = 1, size(m)
            logp = logp + log_choose(m(i), x(i))
         end do
         p = exp(logp)
         return
      end if

      logi = wallenius_log_integral(x, m, odds, precision)
      if (logi /= logi) then
         p = quiet_nan()
         return
      end if
      logp = logi
      do i = 1, size(m)
         logp = logp + log_choose(m(i), x(i))
      end do
      if (logp < log(tiny(1.0_dp))) then
         p = 0.0_dp
      else
         p = exp(logp)
      end if
      p = min(1.0_dp, max(0.0_dp, p))
   end function dmwnchypergeo

   function rmwnchypergeo(nran, m, n, odds, precision) result(draws)
      ! Exact sequential biased-urn simulation.
      integer, intent(in) :: nran, m(:), n
      real(dp), intent(in) :: odds(:)
      real(dp), intent(in), optional :: precision
      integer, allocatable :: draws(:,:)
      integer, allocatable :: remaining(:)
      real(dp), allocatable :: w(:)
      real(dp) :: s, u, c
      integer :: colors, sample, draw, i, selected

      colors = size(m)
      allocate(draws(colors, max(0, nran)))
      if (nran <= 0) return
      if (size(odds) /= colors .or. n < 0 .or. n > sum(m)) then
         draws = -huge(1)
         return
      end if
      allocate(remaining(colors), w(colors))
      do sample = 1, nran
         remaining = m
         draws(:, sample) = 0
         do draw = 1, n
            w = odds * real(remaining, dp)
            s = sum(w)
            if (s <= 0.0_dp) then
               draws(:, sample) = -huge(1)
               exit
            end if
            u = rand_uniform() * s
            c = 0.0_dp
            selected = colors
            do i = 1, colors
               c = c + w(i)
               if (u <= c) then
                  selected = i
                  exit
               end if
            end do
            draws(selected, sample) = draws(selected, sample) + 1
            remaining(selected) = remaining(selected) - 1
         end do
      end do
   end function rmwnchypergeo

   subroutine momentsmwnchypergeo(m, n, odds, mean, variance, precision, combinations)
      ! Exact finite-support summation. Like the upstream moments routine,
      ! cost grows with the number of feasible compositions.
      integer, intent(in) :: m(:), n
      real(dp), intent(in) :: odds(:)
      real(dp), intent(out) :: mean(:), variance(:)
      real(dp), intent(in), optional :: precision
      integer, intent(out), optional :: combinations
      integer, allocatable :: x(:), remaining_capacity(:)
      real(dp), allocatable :: sx(:), sx2(:)
      real(dp) :: sp
      integer :: count, colors, i

      colors = size(m)
      if (size(odds) /= colors .or. size(mean) /= colors .or. size(variance) /= colors) then
         mean = quiet_nan()
         variance = quiet_nan()
         if (present(combinations)) combinations = 0
         return
      end if
      allocate(x(colors), remaining_capacity(colors + 1), sx(colors), sx2(colors))
      remaining_capacity(colors + 1) = 0
      do i = colors, 1, -1
         remaining_capacity(i) = remaining_capacity(i + 1) + m(i)
      end do
      x = 0
      sx = 0.0_dp
      sx2 = 0.0_dp
      sp = 0.0_dp
      count = 0
      call enumerate(1, n)
      if (sp <= 0.0_dp) then
         mean = quiet_nan()
         variance = quiet_nan()
      else
         mean = sx / sp
         variance = max(0.0_dp, sx2 / sp - mean * mean)
      end if
      if (present(combinations)) combinations = count

   contains

      recursive subroutine enumerate(color, left)
         integer, intent(in) :: color, left
         integer :: xi, lo, hi, j
         real(dp) :: pr

         if (color == colors) then
            if (left < 0 .or. left > m(color)) return
            if (odds(color) == 0.0_dp .and. left > 0) return
            x(color) = left
            pr = dmwnchypergeo(x, m, n, odds, precision)
            if (pr /= pr .or. pr <= 0.0_dp) return
            sp = sp + pr
            do j = 1, colors
               sx(j) = sx(j) + real(x(j), dp) * pr
               sx2(j) = sx2(j) + real(x(j) * x(j), dp) * pr
            end do
            if (count < huge(1)) count = count + 1
            return
         end if

         lo = max(0, left - remaining_capacity(color + 1))
         hi = min(m(color), left)
         if (odds(color) == 0.0_dp) hi = min(hi, 0)
         do xi = lo, hi
            x(color) = xi
            call enumerate(color + 1, left - xi)
         end do
      end subroutine enumerate

   end subroutine momentsmwnchypergeo

   function meanmwnchypergeo(m, n, odds, precision) result(mean)
      integer, intent(in) :: m(:), n
      real(dp), intent(in) :: odds(:)
      real(dp), intent(in), optional :: precision
      real(dp) :: mean(size(m)), variance(size(m))
      call momentsmwnchypergeo(m, n, odds, mean, variance, precision)
   end function meanmwnchypergeo

   function varmwnchypergeo(m, n, odds, precision) result(variance)
      integer, intent(in) :: m(:), n
      real(dp), intent(in) :: odds(:)
      real(dp), intent(in), optional :: precision
      real(dp) :: variance(size(m)), mean(size(m))
      call momentsmwnchypergeo(m, n, odds, mean, variance, precision)
   end function varmwnchypergeo

   function oddsmwnchypergeo(mu, m, n, precision) result(odds)
      ! Manly approximation used by BiasedUrn, normalized to a reference color.
      real(dp), intent(in) :: mu(:)
      integer, intent(in) :: m(:), n
      real(dp), intent(in), optional :: precision
      real(dp) :: odds(size(mu))
      integer :: i, c0, ntotal, x1, x2
      real(dp) :: best, d1, d2, a, b

      odds = quiet_nan()
      if (size(mu) /= size(m)) return
      ntotal = sum(m)
      best = 0.0_dp
      c0 = 1
      do i = 1, size(m)
         x1 = max(m(i) + n - ntotal, 0)
         x2 = min(n, m(i))
         d1 = mu(i) - real(x1, dp)
         d2 = real(x2, dp) - mu(i)
         d1 = min(d1, d2)
         if (d1 > best) then
            best = d1
            c0 = i
         end if
      end do
      if (best <= 0.0_dp) return
      odds(c0) = 1.0_dp
      do i = 1, size(m)
         if (i == c0) cycle
         x1 = max(m(i) + n - ntotal, 0)
         x2 = min(n, m(i))
         if (mu(i) < real(x1, dp) .or. mu(i) > real(x2, dp)) then
            odds(i) = quiet_nan()
         else if (mu(i) == real(x1, dp)) then
            odds(i) = 0.0_dp
         else if (mu(i) == real(x2, dp)) then
            odds(i) = positive_inf()
         else
            a = 1.0_dp - mu(i) / real(m(i), dp)
            b = 1.0_dp - mu(c0) / real(m(c0), dp)
            odds(i) = log(a) / log(b)
         end if
      end do
   end function oddsmwnchypergeo

   function nummwnchypergeo(mu_in, n, n_total, odds, precision) result(mout)
      ! Solve Manly's multivariate mean approximation for m.
      real(dp), intent(in) :: mu_in(:), odds(:)
      integer, intent(in) :: n, n_total
      real(dp), intent(in), optional :: precision
      real(dp) :: mout(size(mu_in))
      real(dp) :: mu(size(mu_in)), smu, t, lastt, z, zd, eot, den
      integer :: i, iter
      logical :: bad

      mout = quiet_nan()
      if (size(mu_in) /= size(odds) .or. n <= 0 .or. n > n_total) return
      smu = sum(mu_in)
      if (smu <= 0.0_dp) return
      mu = mu_in * real(n, dp) / smu
      if (n == n_total) then
         mout = mu
         return
      end if
      if (any(odds <= 0.0_dp)) return

      t = -1.0_dp
      do iter = 1, 200
         do
            bad = .false.
            z = 0.0_dp
            zd = 0.0_dp
            do i = 1, size(mu)
               eot = exp(odds(i) * t)
               den = 1.0_dp - eot
               if (den <= 0.0_dp .or. eot <= 0.0_dp) then
                  bad = .true.
                  exit
               end if
               z = z + mu(i) / den
               zd = zd + mu(i) * odds(i) * eot / (den * den)
            end do
            if (.not. bad) exit
            t = 0.125_dp * t
         end do
         lastt = t
         t = t - (z - real(n_total, dp)) / zd
         if (t >= 0.0_dp) t = 0.5_dp * lastt
         if (abs(t - lastt) <= -t * 1.0e-8_dp) exit
      end do
      do i = 1, size(mu)
         mout(i) = mu(i) / (1.0_dp - exp(odds(i) * t))
      end do
   end function nummwnchypergeo

   pure logical function valid_multi(x, m, n, odds) result(ok)
      integer, intent(in) :: x(:), m(:), n
      real(dp), intent(in) :: odds(:)
      integer :: i
      ok = size(x) == size(m) .and. size(m) == size(odds)
      if (.not. ok) return
      ok = n >= 0 .and. n <= sum(m) .and. sum(x) == n .and. all(m >= 0) &
         .and. all(x >= 0) .and. all(odds >= 0.0_dp)
      if (.not. ok) return
      do i = 1, size(m)
         if (x(i) > m(i)) then
            ok = .false.
            return
         end if
      end do
      if (n > sum(pack(m, odds > 0.0_dp))) ok = .false.
   end function valid_multi

end module biasedurn_multiwallenius
