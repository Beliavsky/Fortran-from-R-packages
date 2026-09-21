module ksamples_qn
   use r_kinds, only : dp
   use r_distributions, only : r_qnorm
   use suppdists, only : norm_order
   use ksamples_types, only : qn_result, combined_result, sample_block
   use ksamples_utils, only : score_average_ties, sample_variance
   use ksamples_utils, only : initial_group_labels, next_group_labels, grouped_values
   use ksamples_utils, only : shuffle_real, upper_chisq, multinomial_count, convolve_values
   implicit none
   private
   public :: qn_test, qn_test_combined, qn_raw_statistic

contains

   pure real(dp) function qn_raw_statistic(scores, ns) result(value)
      real(dp), intent(in) :: scores(:) !! Concatenated rank scores grouped according to ns.
      integer, intent(in) :: ns(:) !! Positive sample sizes; sum(ns)=size(scores).
      integer :: i, first, last
      real(dp) :: group_sum
      if (sum(ns) /= size(scores)) error stop 'qn_raw_statistic: invalid group sizes'
      value = 0.0_dp
      first = 1
      do i = 1, size(ns)
         last = first + ns(i) - 1
         group_sum = sum(scores(first:last))
         value = value + group_sum*group_sum/real(ns(i), dp)
         first = last + 1
      end do
      value = anint(1.0e8_dp*value)/1.0e8_dp
   end function qn_raw_statistic

   subroutine qn_test(x, ns, result, score, method, nsim, return_dist)
      real(dp), intent(in) :: x(:) !! Concatenated observations ordered by sample according to ns.
      integer, intent(in) :: ns(:) !! Positive sample sizes; sum(ns)=size(x).
      type(qn_result), intent(out) :: result !! Rank-score statistic, p-values, and optional null distribution.
      character(len=*), intent(in), optional :: score !! "KW", "vdW", or "NS"; default "KW".
      character(len=*), intent(in), optional :: method !! "asymptotic", "simulated", or "exact"; default asymptotic.
      integer, intent(in), optional :: nsim !! Simulation budget and exact-enumeration ceiling; default 10000.
      logical, intent(in), optional :: return_dist !! If true, retain the normalized null statistics.
      real(dp), allocatable :: base_scores(:), averaged(:), work(:), dist(:), ns_scores(:), grouped(:)
      integer, allocatable :: labels(:)
      real(dp) :: raw_obs, obs, smean, svar, ncomb, stat
      integer :: i, n, simulations, hits, visited
      logical :: keep, more
      character(len=10) :: mode
      character(len=3) :: score_mode

      if (size(ns) < 2 .or. any(ns <= 0) .or. sum(ns) /= size(x)) error stop 'qn_test: invalid samples'
      n = size(x)
      score_mode = 'KW'
      if (present(score)) score_mode = adjustl(score)
      allocate(base_scores(n))
      select case (trim(score_mode))
      case ('KW', 'kw')
         do i = 1, n
            base_scores(i) = real(i, dp)
         end do
         result%score = 'KW'
      case ('vdW', 'VDW', 'vdw')
         do i = 1, n
            base_scores(i) = r_qnorm(real(i, dp)/real(n + 1, dp))
         end do
         result%score = 'vdW'
      case ('NS', 'ns')
         call norm_order(n, ns_scores)
         base_scores = ns_scores
         result%score = 'NS'
      case default
         error stop 'qn_test: score must be KW, vdW, or NS'
      end select
      call score_average_ties(x, base_scores, averaged)
      smean = sum(averaged)/real(n, dp)
      svar = sample_variance(averaged)
      if (svar <= 0.0_dp) error stop 'qn_test: zero score variance'
      raw_obs = qn_raw_statistic(averaged, ns)
      obs = (raw_obs - real(n, dp)*smean*smean)/svar
      result%statistic = obs
      result%asymptotic_p = upper_chisq(obs, real(size(ns) - 1, dp))

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
         if (ncomb > real(huge(1), dp)) error stop 'qn_test: exact distribution too large'
         if (keep) allocate(dist(nint(ncomb)))
         call initial_group_labels(ns, labels)
         visited = 0
         more = .true.
         do while (more)
            call grouped_values(averaged, labels, size(ns), grouped)
            visited = visited + 1
            stat = (qn_raw_statistic(grouped, ns) - real(n, dp)*smean*smean)/svar
            if (stat >= obs) hits = hits + 1
            if (keep) dist(visited) = anint(1.0e8_dp*stat)/1.0e8_dp
            call next_group_labels(labels, more)
         end do
         result%randomization_p = real(hits, dp)/real(visited, dp)
         if (keep) result%null_dist = dist(:visited)
      else
         if (keep) allocate(dist(simulations))
         allocate(work(n))
         do visited = 1, simulations
            work = averaged
            call shuffle_real(work)
            call consume_sim(work, visited)
         end do
         result%randomization_p = real(hits, dp)/real(simulations, dp)
         if (keep) result%null_dist = dist
      end if

   contains

      subroutine consume_sim(values, index)
         real(dp), intent(in) :: values(:) !! One random permutation of the pooled scores.
         integer, intent(in) :: index !! Simulation index for storing the null statistic.
         real(dp) :: stat
         stat = (qn_raw_statistic(values, ns) - real(n, dp)*smean*smean)/svar
         if (stat >= obs) hits = hits + 1
         if (keep) dist(index) = anint(1.0e8_dp*stat)/1.0e8_dp
      end subroutine consume_sim
   end subroutine qn_test

   subroutine qn_test_combined(blocks, result, score, method, nsim, return_dist)
      type(sample_block), intent(in) :: blocks(:) !! Independent blocks, each containing concatenated x and sample sizes ns.
      type(combined_result), intent(out) :: result !! Sum of block QN statistics and combined p-values/distribution.
      character(len=*), intent(in), optional :: score !! "KW", "vdW", or "NS"; default "KW".
      character(len=*), intent(in), optional :: method !! "asymptotic", "simulated", or "exact"; default asymptotic.
      integer, intent(in), optional :: nsim !! Simulation budget and joint exact-enumeration ceiling; default 10000.
      logical, intent(in), optional :: return_dist !! If true, retain the combined randomization distribution.
      type(qn_result) :: one
      real(dp), allocatable :: combined(:), tmp(:)
      real(dp) :: joint_count
      integer :: i, simulations, df
      logical :: keep
      character(len=10) :: mode
      character(len=3) :: score_mode

      if (size(blocks) < 1) error stop 'qn_test_combined: no blocks'
      mode = 'asymptotic'
      if (present(method)) mode = adjustl(method)
      score_mode = 'KW'
      if (present(score)) score_mode = adjustl(score)
      simulations = 10000
      if (present(nsim)) simulations = max(1, nsim)
      keep = .false.
      if (present(return_dist)) keep = return_dist
      joint_count = 1.0_dp
      do i = 1, size(blocks)
         if (.not. allocated(blocks(i)%x) .or. .not. allocated(blocks(i)%ns)) error stop 'qn_test_combined: incomplete block'
         joint_count = joint_count*multinomial_count(blocks(i)%ns)
      end do
      if (trim(mode) == 'exact' .and. joint_count > real(simulations, dp)) mode = 'simulated'
      if (trim(mode) /= 'exact' .and. trim(mode) /= 'simulated') mode = 'asymptotic'
      result%method = mode
      result%statistic = 0.0_dp
      df = 0

      do i = 1, size(blocks)
         if (trim(mode) == 'asymptotic') then
            call qn_test(blocks(i)%x, blocks(i)%ns, one, score=score_mode, method='asymptotic')
         else
            call qn_test(blocks(i)%x, blocks(i)%ns, one, score=score_mode, method=mode, &
               nsim=simulations, return_dist=.true.)
         end if
         result%statistic = result%statistic + one%statistic
         df = df + size(blocks(i)%ns) - 1
         if (trim(mode) /= 'asymptotic') then
            if (i == 1) then
               combined = one%null_dist
            else if (trim(mode) == 'simulated') then
               if (size(combined) /= size(one%null_dist)) error stop 'qn_test_combined: simulation lengths differ'
               combined = combined + one%null_dist
            else
               call convolve_values(combined, one%null_dist, tmp)
               call move_alloc(tmp, combined)
            end if
         end if
      end do
      result%asymptotic_p = upper_chisq(result%statistic, real(df, dp))
      if (trim(mode) /= 'asymptotic') then
         result%randomization_p = real(count(combined >= result%statistic), dp)/real(size(combined), dp)
         if (keep) result%null_dist = combined
      end if
   end subroutine qn_test_combined

end module ksamples_qn
