! BiasedUrn-fortran
! Multivariate Fisher's noncentral hypergeometric distribution.
! Upstream BiasedUrn copyright (c) 2002-2024 Agner Fog.
! License: GNU General Public License version 3.
module biasedurn_multifisher
   use biasedurn_kinds, only : dp
   use biasedurn_math, only : log_choose, log_add, log_zero, quiet_nan, &
      positive_inf, sample_log_weights
   implicit none
   private

   public :: dmfnchypergeo, rmfnchypergeo, momentsmfnchypergeo
   public :: meanmfnchypergeo, varmfnchypergeo
   public :: oddsmfnchypergeo, nummfnchypergeo
   public :: minmhypergeo, maxmhypergeo

contains

   function minmhypergeo(m, n) result(xmin)
      integer, intent(in) :: m(:), n
      integer :: xmin(size(m))
      integer :: i, ntotal
      ntotal = sum(m)
      do i = 1, size(m)
         xmin(i) = max(n - ntotal + m(i), 0)
      end do
   end function minmhypergeo

   function maxmhypergeo(m, n) result(xmax)
      integer, intent(in) :: m(:), n
      integer :: xmax(size(m))
      integer :: i
      do i = 1, size(m)
         xmax(i) = min(m(i), n)
      end do
   end function maxmhypergeo

   function mfisher_log_normalizer(m, n, odds) result(logz)
      integer, intent(in) :: m(:), n
      real(dp), intent(in) :: odds(:)
      real(dp) :: logz
      real(dp), allocatable :: dpv(:), newv(:)
      real(dp) :: term
      integer :: i, k, x, kmax

      allocate(dpv(0:n), newv(0:n))
      dpv = log_zero
      dpv(0) = 0.0_dp
      do i = 1, size(m)
         newv = log_zero
         do k = 0, n
            if (dpv(k) <= log_zero / 2.0_dp) cycle
            kmax = min(m(i), n - k)
            do x = 0, kmax
               if (odds(i) == 0.0_dp .and. x > 0) cycle
               term = log_choose(m(i), x)
               if (x > 0) term = term + real(x, dp) * log(odds(i))
               newv(k + x) = log_add(newv(k + x), dpv(k) + term)
            end do
         end do
         dpv = newv
      end do
      logz = dpv(n)
   end function mfisher_log_normalizer

   function dmfnchypergeo(x, m, n, odds, precision) result(p)
      integer, intent(in) :: x(:), m(:), n
      real(dp), intent(in) :: odds(:)
      real(dp), intent(in), optional :: precision
      real(dp) :: p, lognum, logz
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
      lognum = 0.0_dp
      do i = 1, size(m)
         if (x(i) > 0 .and. odds(i) == 0.0_dp) return
         lognum = lognum + log_choose(m(i), x(i))
         if (x(i) > 0) lognum = lognum + real(x(i), dp) * log(odds(i))
      end do
      logz = mfisher_log_normalizer(m, n, odds)
      if (logz <= log_zero / 2.0_dp) then
         p = quiet_nan()
      else
         p = exp(lognum - logz)
      end if
   end function dmfnchypergeo

   function rmfnchypergeo(nran, m, n, odds, precision) result(draws)
      integer, intent(in) :: nran, m(:), n
      real(dp), intent(in) :: odds(:)
      real(dp), intent(in), optional :: precision
      integer, allocatable :: draws(:,:)
      real(dp), allocatable :: suffix(:,:), logw(:)
      real(dp) :: term
      integer :: colors, i, j, r, k, x, xmax, pick

      colors = size(m)
      allocate(draws(colors, max(0, nran)))
      if (nran <= 0) return
      if (size(odds) /= colors) then
         draws = -huge(1)
         return
      end if

      ! suffix(i,r) is log coefficient for colors i..colors selecting r items.
      allocate(suffix(colors + 1, 0:n))
      suffix = log_zero
      suffix(colors + 1, 0) = 0.0_dp
      do i = colors, 1, -1
         do r = 0, n
            do x = 0, min(m(i), r)
               if (odds(i) == 0.0_dp .and. x > 0) cycle
               if (suffix(i + 1, r - x) <= log_zero / 2.0_dp) cycle
               term = log_choose(m(i), x)
               if (x > 0) term = term + real(x, dp) * log(odds(i))
               suffix(i, r) = log_add(suffix(i, r), term + suffix(i + 1, r - x))
            end do
         end do
      end do

      do j = 1, nran
         r = n
         draws(:, j) = 0
         do i = 1, colors - 1
            xmax = min(m(i), r)
            allocate(logw(0:xmax))
            logw = log_zero
            do x = 0, xmax
               if (odds(i) == 0.0_dp .and. x > 0) cycle
               if (suffix(i + 1, r - x) <= log_zero / 2.0_dp) cycle
               term = log_choose(m(i), x)
               if (x > 0) term = term + real(x, dp) * log(odds(i))
               logw(x) = term + suffix(i + 1, r - x)
            end do
            pick = sample_log_weights(logw) - 1
            if (pick < 0) then
               draws(:, j) = -huge(1)
               exit
            end if
            draws(i, j) = pick
            r = r - pick
            deallocate(logw)
         end do
         if (draws(1, j) /= -huge(1)) draws(colors, j) = r
      end do
   end function rmfnchypergeo

   subroutine momentsmfnchypergeo(m, n, odds, mean, variance, precision, combinations)
      integer, intent(in) :: m(:), n
      real(dp), intent(in) :: odds(:)
      real(dp), intent(out) :: mean(:), variance(:)
      real(dp), intent(in), optional :: precision
      integer, intent(out), optional :: combinations
      real(dp), allocatable :: suffix(:,:), prefix(:,:), logw(:)
      real(dp) :: logz, lp, p, s1, s2, term
      integer :: colors, i, j, r, x, count

      colors = size(m)
      if (size(mean) /= colors .or. size(variance) /= colors .or. size(odds) /= colors) then
         mean = quiet_nan()
         variance = quiet_nan()
         if (present(combinations)) combinations = 0
         return
      end if
      logz = mfisher_log_normalizer(m, n, odds)
      mean = 0.0_dp
      variance = 0.0_dp
      if (logz <= log_zero / 2.0_dp) then
         mean = quiet_nan()
         variance = quiet_nan()
         if (present(combinations)) combinations = 0
         return
      end if

      ! For each color build DP for all other colors and obtain its marginal.
      do i = 1, colors
         allocate(logw(0:min(m(i), n)))
         logw = log_zero
         do x = 0, min(m(i), n)
            if (odds(i) == 0.0_dp .and. x > 0) cycle
            lp = other_log_normalizer(m, n - x, odds, i)
            if (lp <= log_zero / 2.0_dp) cycle
            term = log_choose(m(i), x)
            if (x > 0) term = term + real(x, dp) * log(odds(i))
            logw(x) = term + lp - logz
         end do
         s1 = 0.0_dp
         s2 = 0.0_dp
         do x = 0, ubound(logw, 1)
            if (logw(x) <= log_zero / 2.0_dp) cycle
            p = exp(logw(x))
            s1 = s1 + real(x, dp) * p
            s2 = s2 + real(x * x, dp) * p
         end do
         mean(i) = s1
         variance(i) = max(0.0_dp, s2 - s1 * s1)
         deallocate(logw)
      end do

      if (present(combinations)) then
         ! Count feasible compositions by an integer DP.
         count = count_combinations(m, n)
         combinations = count
      end if
   end subroutine momentsmfnchypergeo

   function meanmfnchypergeo(m, n, odds, precision) result(mean)
      integer, intent(in) :: m(:), n
      real(dp), intent(in) :: odds(:)
      real(dp), intent(in), optional :: precision
      real(dp) :: mean(size(m)), variance(size(m))
      call momentsmfnchypergeo(m, n, odds, mean, variance, precision)
   end function meanmfnchypergeo

   function varmfnchypergeo(m, n, odds, precision) result(variance)
      integer, intent(in) :: m(:), n
      real(dp), intent(in) :: odds(:)
      real(dp), intent(in), optional :: precision
      real(dp) :: variance(size(m)), mean(size(m))
      call momentsmfnchypergeo(m, n, odds, mean, variance, precision)
   end function varmfnchypergeo

   function oddsmfnchypergeo(mu, m, n, precision) result(odds)
      ! Cornfield approximation used by BiasedUrn, normalized to one
      ! reference color chosen to be furthest from its support boundary.
      real(dp), intent(in) :: mu(:)
      integer, intent(in) :: m(:), n
      real(dp), intent(in), optional :: precision
      real(dp) :: odds(size(mu))
      integer :: i, c0, ntotal, x1, x2
      real(dp) :: best, d1, d2

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
            odds(i) = mu(i) * (real(m(c0), dp) - mu(c0)) / &
               (mu(c0) * (real(m(i), dp) - mu(i)))
         end if
      end do
   end function oddsmfnchypergeo

   function nummfnchypergeo(mu_in, n, n_total, odds, precision) result(mout)
      real(dp), intent(in) :: mu_in(:), odds(:)
      integer, intent(in) :: n, n_total
      real(dp), intent(in), optional :: precision
      real(dp) :: mout(size(mu_in))
      real(dp) :: mu(size(mu_in)), smu, r, z, zd, last
      integer :: i, iter

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

      r = 1.0_dp
      do iter = 1, 200
         last = r
         z = 0.0_dp
         zd = 0.0_dp
         do i = 1, size(mu)
            z = z + mu(i) * (1.0_dp + 1.0_dp / (r * odds(i)))
            zd = zd - mu(i) / (odds(i) * r * r)
         end do
         r = r - (z - real(n_total, dp)) / zd
         if (r <= 0.0_dp) r = 0.5_dp * last
         if (abs(r - last) <= r * 1.0e-8_dp) exit
      end do
      do i = 1, size(mu)
         mout(i) = mu(i) * (r * odds(i) + 1.0_dp) / (r * odds(i))
      end do
   end function nummfnchypergeo

   function other_log_normalizer(m, n, odds, skip) result(logz)
      integer, intent(in) :: m(:), n, skip
      real(dp), intent(in) :: odds(:)
      real(dp) :: logz
      real(dp), allocatable :: dpv(:), newv(:)
      real(dp) :: term
      integer :: i, k, x
      if (n < 0) then
         logz = log_zero
         return
      end if
      allocate(dpv(0:n), newv(0:n))
      dpv = log_zero
      dpv(0) = 0.0_dp
      do i = 1, size(m)
         if (i == skip) cycle
         newv = log_zero
         do k = 0, n
            if (dpv(k) <= log_zero / 2.0_dp) cycle
            do x = 0, min(m(i), n - k)
               if (odds(i) == 0.0_dp .and. x > 0) cycle
               term = log_choose(m(i), x)
               if (x > 0) term = term + real(x, dp) * log(odds(i))
               newv(k + x) = log_add(newv(k + x), dpv(k) + term)
            end do
         end do
         dpv = newv
      end do
      logz = dpv(n)
   end function other_log_normalizer

   integer function count_combinations(m, n) result(count)
      integer, intent(in) :: m(:), n
      integer, allocatable :: dpv(:), newv(:)
      integer :: i, k, x
      allocate(dpv(0:n), newv(0:n))
      dpv = 0
      dpv(0) = 1
      do i = 1, size(m)
         newv = 0
         do k = 0, n
            if (dpv(k) == 0) cycle
            do x = 0, min(m(i), n - k)
               if (newv(k + x) > huge(1) - dpv(k)) then
                  newv(k + x) = huge(1)
               else
                  newv(k + x) = newv(k + x) + dpv(k)
               end if
            end do
         end do
         dpv = newv
      end do
      count = dpv(n)
   end function count_combinations

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

end module biasedurn_multifisher
