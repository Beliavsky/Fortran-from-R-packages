module ksamples_contingency
   use r_kinds, only : dp
   use ksamples_types, only : contingency_block, contingency_result, combined_result
   use ksamples_utils, only : combination_count, log_choose, upper_chisq
   use ksamples_utils, only : shuffle_integer, convolve_distribution, relative_frequencies
   implicit none
   private
   public :: contingency2xt, contingency2xt_comb

contains

   pure real(dp) function kw_from_counts(avec, bvec) result(kw)
      integer, intent(in) :: avec(:) !! Counts in response row A for each treatment column.
      integer, intent(in) :: bvec(:) !! Counts in response row B for each treatment column.
      integer :: m, n, total
      integer, allocatable :: d(:)
      real(dp) :: delta
      if (size(avec) /= size(bvec) .or. size(avec) < 2) error stop 'kw_from_counts: invalid table'
      if (any(avec < 0) .or. any(bvec < 0)) error stop 'kw_from_counts: negative count'
      d = avec + bvec
      if (any(d <= 0)) error stop 'kw_from_counts: empty column'
      m = sum(avec)
      n = sum(bvec)
      total = m + n
      if (m <= 0 .or. n <= 0) error stop 'kw_from_counts: empty row'
      delta = sum(real(avec*avec, dp)/real(d, dp))
      kw = real(total*(total - 1), dp)/real(m*n, dp)*(delta - real(m*m, dp)/real(total, dp))
   end function kw_from_counts

   subroutine contingency2xt(avec, bvec, result, method, nsim, return_dist)
      integer, intent(in) :: avec(:) !! Counts in response row A, one per treatment column.
      integer, intent(in) :: bvec(:) !! Counts in response row B, conforming with avec.
      type(contingency_result), intent(out) :: result !! KW statistic, p-values, and optional tabulated null distribution.
      character(len=*), intent(in), optional :: method !! "asymptotic", "simulated", or "exact"; default asymptotic.
      integer, intent(in), optional :: nsim !! Simulation budget and exact-composition ceiling; default 1000000.
      logical, intent(in), optional :: return_dist !! If true, retain a tabulated exact or empirical null distribution.
      integer, allocatable :: d(:), binary(:), sim_a(:)
      real(dp), allocatable :: vals(:), probs(:), simvals(:), support(:), probability(:)
      real(dp) :: observed, delta_obs, ncomb, denom_log, exact_p
      integer :: m, n, total, tnum, simulations, count_valid, isim, offset, j
      logical :: keep
      character(len=10) :: mode

      if (size(avec) /= size(bvec) .or. size(avec) < 2) error stop 'contingency2xt: nonconforming table'
      if (any(avec < 0) .or. any(bvec < 0)) error stop 'contingency2xt: negative count'
      d = avec + bvec
      if (any(d <= 0)) error stop 'contingency2xt: empty column'
      m = sum(avec)
      n = sum(bvec)
      total = m + n
      tnum = size(avec)
      if (m <= 0 .or. n <= 0) error stop 'contingency2xt: empty response row'
      observed = kw_from_counts(avec, bvec)
      result%statistic = anint(1.0e6_dp*observed)/1.0e6_dp
      result%asymptotic_p = upper_chisq(observed, real(tnum - 1, dp))
      mode = 'asymptotic'
      if (present(method)) mode = adjustl(method)
      simulations = 1000000
      if (present(nsim)) simulations = max(100, nsim)
      keep = .false.
      if (present(return_dist)) keep = return_dist
      ncomb = combination_count(m + tnum - 1, tnum - 1)
      if (trim(mode) == 'exact' .and. ncomb > real(simulations, dp)) mode = 'simulated'
      if (trim(mode) /= 'exact' .and. trim(mode) /= 'simulated') mode = 'asymptotic'
      result%method = mode
      if (trim(mode) == 'asymptotic') return

      delta_obs = sum(real(avec*avec, dp)/real(d, dp))
      if (trim(mode) == 'exact') then
         allocate(vals(nint(ncomb)), probs(nint(ncomb)))
         vals = 0.0_dp
         probs = 0.0_dp
         count_valid = 0
         exact_p = 0.0_dp
         denom_log = log_choose(total, m)
         allocate(sim_a(tnum))
         call enumerate_composition(1, m)
         result%randomization_p = exact_p
         if (keep) call compress_weighted(vals(:count_valid), probs(:count_valid), result%support, result%probability)
      else
         allocate(binary(total), sim_a(tnum), simvals(simulations))
         binary = 0
         binary(:m) = 1
         do isim = 1, simulations
            call shuffle_integer(binary)
            offset = 0
            do j = 1, tnum
               sim_a(j) = sum(binary(offset + 1:offset + d(j)))
               offset = offset + d(j)
            end do
            simvals(isim) = kw_from_counts(sim_a, d - sim_a)
         end do
         result%randomization_p = real(count(simvals >= observed), dp)/real(simulations, dp)
         if (keep) then
            call relative_frequencies(simvals, support, probability, digits=6)
            result%support = support
            result%probability = probability
         end if
      end if

   contains
      recursive subroutine enumerate_composition(col, remaining)
         integer, intent(in) :: col !! Current table column being assigned a row-A count.
         integer, intent(in) :: remaining !! Row-A count remaining for columns col:tnum.
         integer :: a, i
         real(dp) :: delta, prob, kw
         if (col == tnum) then
            sim_a(col) = remaining
            if (remaining < 0 .or. remaining > d(col)) return
            delta = 0.0_dp
            prob = 0.0_dp
            do i = 1, tnum
               delta = delta + real(sim_a(i)*sim_a(i), dp)/real(d(i), dp)
               prob = prob + log_choose(d(i), sim_a(i))
            end do
            prob = exp(prob - denom_log)
            if (delta >= delta_obs) exact_p = exact_p + prob
            if (keep) then
               count_valid = count_valid + 1
               kw = real(total*(total - 1), dp)/real(m*n, dp)*(delta - real(m*m, dp)/real(total, dp))
               vals(count_valid) = anint(1.0e6_dp*kw)/1.0e6_dp
               probs(count_valid) = prob
            end if
            return
         end if
         do a = 0, remaining
            if (a <= d(col)) then
               sim_a(col) = a
               call enumerate_composition(col + 1, remaining - a)
            end if
         end do
      end subroutine enumerate_composition
   end subroutine contingency2xt

   subroutine contingency2xt_comb(blocks, result, method, nsim, return_dist)
      type(contingency_block), intent(in) :: blocks(:) !! Independent 2 x t blocks; t may vary by block.
      type(combined_result), intent(out) :: result !! Sum of block KW statistics and combined p-values/distribution.
      character(len=*), intent(in), optional :: method !! "asymptotic", "simulated", or "exact"; default asymptotic.
      integer, intent(in), optional :: nsim !! Simulation budget and joint exact-enumeration ceiling; default 10000.
      logical, intent(in), optional :: return_dist !! If true, retain a tabulated combined null distribution.
      type(contingency_result) :: one
      real(dp), allocatable :: sx(:), sp(:), tx(:), tp(:), simsum(:), blocksim(:)
      integer, allocatable :: d(:), binary(:), sim_a(:)
      real(dp) :: joint_count
      integer :: i, j, isim, m, n, total, offset, simulations, df
      logical :: keep
      character(len=10) :: mode

      if (size(blocks) < 2) error stop 'contingency2xt_comb: need at least two blocks'
      mode = 'asymptotic'
      if (present(method)) mode = adjustl(method)
      simulations = 10000
      if (present(nsim)) simulations = max(100, nsim)
      keep = .false.
      if (present(return_dist)) keep = return_dist
      joint_count = 1.0_dp
      do i = 1, size(blocks)
         if (.not. allocated(blocks(i)%avec) .or. .not. allocated(blocks(i)%bvec)) &
            error stop 'contingency2xt_comb: incomplete block'
         m = sum(blocks(i)%avec)
         joint_count = joint_count*combination_count(m + size(blocks(i)%avec) - 1, size(blocks(i)%avec) - 1)
      end do
      if (trim(mode) == 'exact' .and. joint_count > real(simulations, dp)) mode = 'simulated'
      if (trim(mode) /= 'exact' .and. trim(mode) /= 'simulated') mode = 'asymptotic'
      result%method = mode
      result%statistic = 0.0_dp
      df = 0

      if (trim(mode) == 'simulated') then
         allocate(simsum(simulations))
         simsum = 0.0_dp
      end if
      do i = 1, size(blocks)
         df = df + size(blocks(i)%avec) - 1
         if (trim(mode) == 'asymptotic') then
            call contingency2xt(blocks(i)%avec, blocks(i)%bvec, one, method='asymptotic')
            result%statistic = result%statistic + one%statistic
         else if (trim(mode) == 'exact') then
            call contingency2xt(blocks(i)%avec, blocks(i)%bvec, one, method='exact', nsim=simulations, return_dist=.true.)
            result%statistic = result%statistic + one%statistic
            if (i == 1) then
               sx = one%support
               sp = one%probability
            else
               call convolve_distribution(sx, sp, one%support, one%probability, tx, tp)
               call move_alloc(tx, sx)
               call move_alloc(tp, sp)
            end if
         else
            result%statistic = result%statistic + kw_from_counts(blocks(i)%avec, blocks(i)%bvec)
            d = blocks(i)%avec + blocks(i)%bvec
            m = sum(blocks(i)%avec)
            n = sum(blocks(i)%bvec)
            total = m + n
            allocate(binary(total), sim_a(size(d)), blocksim(simulations))
            binary = 0
            binary(:m) = 1
            do isim = 1, simulations
               call shuffle_integer(binary)
               offset = 0
               do j = 1, size(d)
                  sim_a(j) = sum(binary(offset + 1:offset + d(j)))
                  offset = offset + d(j)
               end do
               blocksim(isim) = kw_from_counts(sim_a, d - sim_a)
            end do
            simsum = simsum + blocksim
            deallocate(binary, sim_a, blocksim, d)
         end if
      end do
      result%asymptotic_p = upper_chisq(result%statistic, real(df, dp))
      if (trim(mode) == 'exact') then
         result%randomization_p = sum(sp, mask=sx >= result%statistic)
         if (keep) then
            result%support = sx
            result%probability = sp
         end if
      else if (trim(mode) == 'simulated') then
         result%randomization_p = real(count(simsum >= result%statistic), dp)/real(simulations, dp)
         if (keep) call relative_frequencies(simsum, result%support, result%probability, digits=6)
      end if
   end subroutine contingency2xt_comb

   pure subroutine compress_weighted(values, weights, support, probability)
      real(dp), intent(in) :: values(:) !! Rounded support candidates.
      real(dp), intent(in) :: weights(:) !! Nonnegative weights corresponding to values.
      real(dp), allocatable, intent(out) :: support(:) !! Sorted unique support.
      real(dp), allocatable, intent(out) :: probability(:) !! Sum of weights at each support point.
      real(dp), allocatable :: sv(:), sw(:)
      real(dp) :: v, w
      integer :: i, n, pos
      if (size(values) /= size(weights)) error stop 'compress_weighted: shape mismatch'
      allocate(sv(size(values)), sw(size(values)))
      n = 0
      do i = 1, size(values)
         v = values(i)
         w = weights(i)
         pos = 1
         do while (pos <= n)
            if (sv(pos) >= v) exit
            pos = pos + 1
         end do
         if (pos <= n) then
            if (sv(pos) == v) then
               sw(pos) = sw(pos) + w
               cycle
            end if
         end if
            if (pos <= n) then
               sv(pos + 1:n + 1) = sv(pos:n)
               sw(pos + 1:n + 1) = sw(pos:n)
            end if
            sv(pos) = v
            sw(pos) = w
            n = n + 1
      end do
      allocate(support(n), probability(n))
      support = sv(:n)
      probability = sw(:n)
   end subroutine compress_weighted

end module ksamples_contingency
