! SPDX-License-Identifier: GPL-3.0-or-later
program test_reference
   use monobin
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   character(len=160) :: name
   character(len=3) :: method
   real(dp), allocatable :: x(:), y(:), specials(:), pred(:), expected(:,:), w(:)
   integer, allocatable :: membership(:), expected_member(:), counts(:)
   type(bin_options) :: options
   type(bin_result) :: model
   integer :: index_unit, input, ios, n, ns, separate, trend, nb, nsc, i, j, failures, tested
   failures = 0
   tested = 0
   open(newunit=index_unit, file='test/fixtures/index.txt', status='old', action='read')
   do
      read(index_unit, '(a)', iostat=ios) name
      if (ios /= 0) exit
      open(newunit=input, file='test/fixtures/' // trim(name), status='old', action='read')
      options = bin_options()
      read(input, *) n, method, options%target, options%direction, options%g, separate, trend, &
         options%min_pct_obs, options%min_avg_rate, options%woe_gap, options%p_value, ns
      options%separate_special = separate == 1
      options%woe_trend = trend == 1
      allocate(x(n), y(n), specials(ns), expected_member(n))
      if (ns > 0) read(input, *) specials
      do i = 1, n
         read(input, *) x(i), y(i)
      end do
      read(input, *) nb, nsc
      allocate(counts(nb), expected(8, nb))
      do i = 1, nb
         read(input, *) counts(i), expected(:, i)
      end do
      read(input, *) expected_member
      close(input)
      call fit_bins(x, y, method, model, options, specials)
      tested = tested + 1
      if (model%status /= bin_ok .or. size(model%bins) /= nb .or. model%n_special /= nsc) then
         print *, 'STRUCTURE FAIL ', trim(name), model%status, size(model%bins), nb
         failures = failures + 1
      else
         do i = 1, nb
            if (model%bins(i)%count /= counts(i)) call fail('count')
            call near(model%bins(i)%y_sum, expected(1, i), 'y sum')
            call near(model%bins(i)%y_mean, expected(2, i), 'y mean')
            ! Special x summaries depend on order-sensitive +/-Inf and NaN reductions.
            if (i > nsc) then
               call near(model%bins(i)%x_mean, expected(3, i), 'x mean')
               call near(model%bins(i)%x_min, expected(4, i), 'x min')
               call near(model%bins(i)%x_max, expected(5, i), 'x max')
            end if
            call near(model%bins(i)%woe, expected(6, i), 'woe')
            call near(model%bins(i)%iv, expected(7, i), 'iv')
            if (method == 'sts') call near(model%bins(i)%p_value, expected(8, i), 'adjacent p value')
         end do
         if (any(model%membership /= expected_member)) call fail('membership')
      end if
      if (model%status == bin_ok) then
         ! Independent calibration, assignment, and normalized-share checks.
         if (sum(model%bins%count) /= count(.not. ieee_is_nan(y))) call fail('retained counts')
         call near(sum(model%bins%y_sum), sum(pack(y, .not. ieee_is_nan(y))), 'response conservation')
         call apply_bins(model, x, membership, pred, w)
         do i = 1, size(model%bins)
            if (count(model%membership == i) /= model%bins(i)%count) call fail('bin membership count')
            call near(sum(pack(y, model%membership == i)), model%bins(i)%y_sum, 'bin calibration')
         end do
         do i = 1, n
            if (ieee_is_nan(y(i))) cycle
            if (membership(i) /= model%membership(i)) call fail('apply consistency')
            call near(pred(i), model%bins(membership(i))%y_mean, 'prediction calibration')
         end do
         if (all(ieee_is_finite(model%bins%dist_good))) call near(sum(model%bins%dist_good), 1.0_dp, 'share normalization')
         call near(sum(model%bins%dist_bad), 1.0_dp, 'share normalization')
         do i = nsc + 1, size(model%bins) - 1
            if (model%bins(i)%x_max >= model%bins(i + 1)%x_min) call fail('ordered disjoint bins')
         end do
      end if
      deallocate(x, y, specials, counts, expected, expected_member)
   end do
   close(index_unit)
   print *, 'Genuine R fixtures:', tested, 'failures:', failures
   if (failures > 0) error stop 1
contains
   subroutine fail(what)
      character(len=*), intent(in) :: what !! Name of the violated numeric or structural test invariant.
      failures = failures + 1
      if (failures < 60) print *, 'FAIL ', trim(name), ' : ', what
   end subroutine fail
   subroutine near(actual, wanted, what)
      real(dp), intent(in) :: actual !! Actual statistic from the translated fit.
      real(dp), intent(in) :: wanted !! Genuine upstream R statistic or independently conserved total.
      character(len=*), intent(in) :: what !! Human-readable statistic name for a failed comparison.
      if (ieee_is_nan(wanted)) then
         if (.not. ieee_is_nan(actual)) call fail(what)
      else if (.not. ieee_is_finite(wanted)) then
         if (actual /= wanted) call fail(what)
      else if (.not. ieee_is_finite(actual)) then
         call fail(what)
      else if (abs(actual - wanted) > 3.0e-9_dp * max(1.0_dp, abs(wanted))) then
         call fail(what)
         if (failures < 60) print *, actual, wanted
      end if
   end subroutine near
end program test_reference
