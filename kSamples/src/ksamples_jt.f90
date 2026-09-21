module ksamples_jt
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf, ieee_quiet_nan
   use r_kinds, only : dp
   use ksamples_types, only : jt_result
   use ksamples_utils, only : average_ranks, shuffle_real
   use ksamples_utils, only : initial_group_labels, next_group_labels, grouped_values
   use ksamples_utils, only : upper_normal, multinomial_count, combination_count, sort_real_in_place
   implicit none
   private
   public :: jt_test, jt_statistic, harding_distribution, djt, pjt, qjt

contains

   pure real(dp) function jt_statistic(rx, ns) result(value)
      real(dp), intent(in) :: rx(:) !! Concatenated rank scores grouped according to ns.
      integer, intent(in) :: ns(:) !! Positive sample sizes; sum(ns) must equal size(rx).
      integer :: i, j, m, n, m_start, m_end, n_start, n_end
      if (sum(ns) /= size(rx)) error stop 'jt_statistic: group sizes do not sum to data length'
      value = 0.0_dp
      m_start = 1
      do i = 1, size(ns) - 1
         m_end = m_start + ns(i) - 1
         n_start = m_end + 1
         do j = i + 1, size(ns)
            n_end = n_start + ns(j) - 1
            do n = n_start, n_end
               do m = m_start, m_end
                  if (rx(m) < rx(n)) then
                     value = value + 1.0_dp
                  else if (rx(m) == rx(n)) then
                     value = value + 0.5_dp
                  end if
               end do
            end do
            n_start = n_end + 1
         end do
         m_start = m_end + 1
      end do
   end function jt_statistic

   pure subroutine jt_mean_sigma(rx, ns, mean_value, sigma)
      real(dp), intent(in) :: rx(:) !! Average ranks of all observations.
      integer, intent(in) :: ns(:) !! Sample sizes defining the ordered groups.
      real(dp), intent(out) :: mean_value !! Null mean of the JT statistic.
      real(dp), intent(out) :: sigma !! Tie-corrected null standard deviation.
      real(dp), allocatable :: sorted(:)
      integer, allocatable :: ties(:)
      integer :: i, j, nt, n
      real(dp) :: x1, x2, x3, a1, a2, a3, ss2, ss3, tt2, tt3

      n = size(rx)
      sorted = rx
      call sort_real_in_place(sorted)
      allocate(ties(max(1, n)))
      nt = 0
      i = 1
      do while (i <= n)
         j = i
         do while (j < n)
            if (sorted(j + 1) /= sorted(i)) exit
            j = j + 1
         end do
         nt = nt + 1
         ties(nt) = j - i + 1
         i = j + 1
      end do

      x1 = real(n*(n - 1)*(2*n + 5), dp)
      x2 = sum(real(ns*(ns - 1)*(2*ns + 5), dp))
      x3 = sum(real(ties(:nt)*(ties(:nt) - 1)*(2*ties(:nt) + 5), dp))
      a1 = (x1 - x2 - x3)/72.0_dp
      ss3 = sum(real(ns*(ns - 1)*(ns - 2), dp))
      tt3 = sum(real(ties(:nt)*(ties(:nt) - 1)*(ties(:nt) - 2), dp))
      if (n > 2) then
         a2 = ss3*tt3/(36.0_dp*real(n*(n - 1)*(n - 2), dp))
      else
         a2 = 0.0_dp
      end if
      ss2 = sum(real(ns*(ns - 1), dp))
      tt2 = sum(real(ties(:nt)*(ties(:nt) - 1), dp))
      if (n > 1) then
         a3 = ss2*tt2/(8.0_dp*real(n*(n - 1), dp))
      else
         a3 = 0.0_dp
      end if
      sigma = sqrt(max(0.0_dp, a1 + a2 + a3))
      mean_value = 0.0_dp
      do i = 1, size(ns) - 1
         do j = i + 1, size(ns)
            mean_value = mean_value + 0.5_dp*real(ns(i)*ns(j), dp)
         end do
      end do
   end subroutine jt_mean_sigma

   subroutine jt_test(x, ns, result, method, nsim, return_dist)
      real(dp), intent(in) :: x(:) !! Concatenated observations, ordered by sample according to ns.
      integer, intent(in) :: ns(:) !! Positive sample sizes; sum(ns)=size(x).
      type(jt_result), intent(out) :: result !! JT statistic, null moments, p-values, and optional null distribution.
      character(len=*), intent(in), optional :: method !! "asymptotic", "simulated", or "exact"; default asymptotic.
      integer, intent(in), optional :: nsim !! Simulation budget and exact-enumeration ceiling; default 10000.
      logical, intent(in), optional :: return_dist !! If true, retain the simulated or exact null statistics.
      real(dp), allocatable :: ranks(:), work(:), dist(:), grouped(:)
      integer, allocatable :: labels(:)
      real(dp) :: ncomb, obs, stat
      integer :: simulations, hits, visited
      logical :: keep, more
      character(len=10) :: mode

      if (size(ns) < 2 .or. any(ns <= 0) .or. sum(ns) /= size(x)) error stop 'jt_test: invalid samples'
      call average_ranks(x, ranks)
      obs = jt_statistic(ranks, ns)
      call jt_mean_sigma(ranks, ns, result%mean, result%sigma)
      result%statistic = obs
      if (result%sigma > 0.0_dp) then
         result%asymptotic_p = upper_normal((obs - result%mean)/result%sigma)
      else
         result%asymptotic_p = 1.0_dp
      end if

      mode = 'asymptotic'
      if (present(method)) mode = adjustl(method)
      simulations = 10000
      if (present(nsim)) simulations = max(1, nsim)
      keep = .false.
      if (present(return_dist)) keep = return_dist
      ncomb = multinomial_count(ns)
      if (trim(mode) == 'exact' .and. ncomb > real(simulations, dp)) mode = 'simulated'
      if (trim(mode) /= 'exact' .and. trim(mode) /= 'simulated') mode = 'asymptotic'
      result%method = mode
      if (trim(mode) == 'asymptotic') return

      hits = 0
      if (trim(mode) == 'exact') then
         if (ncomb > real(huge(1), dp)) error stop 'jt_test: exact distribution too large'
         if (keep) allocate(dist(nint(ncomb)))
         call initial_group_labels(ns, labels)
         visited = 0
         more = .true.
         do while (more)
            call grouped_values(ranks, labels, size(ns), grouped)
            visited = visited + 1
            stat = jt_statistic(grouped, ns)
            if (stat >= obs) hits = hits + 1
            if (keep) dist(visited) = stat
            call next_group_labels(labels, more)
         end do
         result%randomization_p = real(hits, dp)/real(visited, dp)
         if (keep) result%null_dist = dist(:visited)
      else
         if (keep) allocate(dist(simulations))
         allocate(work(size(ranks)))
         do visited = 1, simulations
            work = ranks
            call shuffle_real(work)
            call consume_sim(work, visited)
         end do
         result%randomization_p = real(hits, dp)/real(simulations, dp)
         if (keep) result%null_dist = dist
      end if

   contains

      subroutine consume_sim(values, index)
         real(dp), intent(in) :: values(:) !! One random permutation of pooled rank scores.
         integer, intent(in) :: index !! Simulation index used when retaining the null distribution.
         real(dp) :: stat
         stat = jt_statistic(values, ns)
         if (stat >= obs) hits = hits + 1
         if (keep) dist(index) = stat
      end subroutine consume_sim
   end subroutine jt_test

   pure subroutine harding_distribution(nn, freq)
      integer, intent(in) :: nn(:) !! Positive sample sizes for the exact Jonckheere-Terpstra distribution.
      real(dp), allocatable, intent(out) :: freq(:) !! Probabilities for statistic values 0:L.
      integer, allocatable :: sizes(:), nvec(:)
      integer :: k, i, j, key, lmax, mid, m, n, pmax, qmax, t, u, s
      real(dp) :: denom
      if (size(nn) < 2 .or. any(nn <= 0)) error stop 'harding_distribution: invalid sample sizes'
      sizes = nn
      do i = 2, size(sizes)
         key = sizes(i)
         j = i - 1
         do while (j >= 1)
            if (sizes(j) <= key) exit
            sizes(j + 1) = sizes(j)
            j = j - 1
         end do
         sizes(j + 1) = key
      end do
      k = size(sizes)
      allocate(nvec(k))
      nvec(k) = sizes(k)
      do i = k - 1, 1, -1
         nvec(i) = sizes(i) + nvec(i + 1)
      end do
      lmax = sum(sizes(:k - 1)*nvec(2:k))
      allocate(freq(lmax + 1))
      freq = 0.0_dp
      freq(1) = 1.0_dp
      mid = lmax/2
      do i = 2, k
         m = nvec(i - 1) - nvec(i)
         n = nvec(i)
         if (n + 1 <= mid) then
            pmax = min(m + n, mid)
            do t = n + 1, pmax
               do u = mid, t, -1
                  freq(u + 1) = freq(u + 1) - freq(u - t + 1)
               end do
            end do
         end if
         qmax = min(mid, m)
         do s = 1, qmax
            do u = s, mid
               freq(u + 1) = freq(u + 1) + freq(u - s + 1)
            end do
         end do
         denom = combination_count(m + n, m)
         freq = freq/denom
      end do
      if (mod(lmax, 2) == 0) then
         do i = 1, mid
            freq(mid + i + 1) = freq(mid - i + 1)
         end do
      else
         do i = 1, mid + 1
            freq(mid + i + 1) = freq(mid + 2 - i)
         end do
      end if
   end subroutine harding_distribution

   pure subroutine djt(x, nn, d)
      real(dp), intent(in) :: x(:) !! JT statistic values at which the exact probability mass is requested.
      integer, intent(in) :: nn(:) !! Sample sizes defining the exact null distribution.
      real(dp), allocatable, intent(out) :: d(:) !! Exact probability mass at each x; zero off the integer support.
      real(dp), allocatable :: freq(:)
      integer :: i, ix
      call harding_distribution(nn, freq)
      allocate(d(size(x)))
      d = 0.0_dp
      do i = 1, size(x)
         ix = nint(x(i))
         if (abs(x(i) - real(ix, dp)) <= 4.0_dp*epsilon(1.0_dp) .and. ix >= 0 .and. ix < size(freq)) then
            d(i) = freq(ix + 1)
         end if
      end do
   end subroutine djt

   pure subroutine pjt(x, nn, p)
      real(dp), intent(in) :: x(:) !! JT values; values are floored as in the R implementation.
      integer, intent(in) :: nn(:) !! Sample sizes defining the exact null distribution.
      real(dp), allocatable, intent(out) :: p(:) !! Exact cumulative probabilities P(JT <= floor(x)).
      real(dp), allocatable :: freq(:), cdf(:)
      integer :: i, ix
      call harding_distribution(nn, freq)
      allocate(cdf(size(freq)), p(size(x)))
      cdf(1) = freq(1)
      do i = 2, size(freq)
         cdf(i) = cdf(i - 1) + freq(i)
      end do
      do i = 1, size(x)
         ix = floor(x(i))
         if (ix < 0) then
            p(i) = 0.0_dp
         else if (ix >= size(freq) - 1) then
            p(i) = 1.0_dp
         else
            p(i) = cdf(ix + 1)
         end if
      end do
   end subroutine pjt

   pure subroutine qjt(prob, nn, quantile)
      real(dp), intent(in) :: prob(:) !! Cumulative probabilities in [0,1].
      integer, intent(in) :: nn(:) !! Sample sizes defining the exact null distribution.
      real(dp), allocatable, intent(out) :: quantile(:) !! Exact inverse-CDF values; negative infinity for p<=0 and NaN above 1.
      real(dp), allocatable :: freq(:)
      real(dp) :: cdf
      integer :: i, j
      call harding_distribution(nn, freq)
      allocate(quantile(size(prob)))
      do i = 1, size(prob)
         if (prob(i) <= 0.0_dp) then
            quantile(i) = ieee_value(0.0_dp, ieee_negative_inf)
         else if (prob(i) > 1.0_dp) then
            quantile(i) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            cdf = 0.0_dp
            quantile(i) = real(size(freq) - 1, dp)
            do j = 1, size(freq)
               cdf = cdf + freq(j)
               if (cdf >= prob(i)) then
                  quantile(i) = real(j - 1, dp)
                  exit
               end if
            end do
         end if
      end do
   end subroutine qjt

end module ksamples_jt
