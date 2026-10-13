! SPDX-License-Identifier: GPL-3.0-or-later
program test_invariants
   use monobin
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, &
      ieee_negative_inf, ieee_is_nan, ieee_is_finite
   implicit none
   character(len=3), parameter :: methods(7) = ['pct', 'iso', 'woe', 'cum', 'sts', 'ndr', 'mdt']
   real(dp) :: x(400), y(400), nan, inf
   real(dp), allocatable :: pred(:), woe(:)
   integer, allocatable :: group(:)
   type(bin_result) :: model
   type(bin_options) :: options
   integer :: i, j, m, direction, failures
   failures = 0
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   do i = 1, size(x)
      x(i) = real((i - 1) / 20, dp)
      y(i) = merge(1.0_dp, 0.0_dp, modulo(13 * i, 100) < 15 + int(3.0_dp * x(i)))
   end do
   do direction = -1, 1, 2
      options = bin_options(direction=direction)
      do m = 1, 7
         call fit_bins(x, y, methods(m), model, options)
         call require(model%status == bin_ok, 'valid finite input')
         do i = model%n_special + 1, size(model%bins) - 1
            call require(real(direction, dp) * (model%bins(i + 1)%y_mean - model%bins(i)%y_mean) &
               >= -1.0e-12_dp, 'forced monotonic response')
         end do
         call apply_bins(model, [-1000.0_dp, 1.0e6_dp], group, pred, woe)
         call require(group(1) == 1 .and. group(2) == size(model%bins), 'tail extrapolation')
         if (size(model%breaks) > 0) then
            call apply_bins(model, model%breaks, group)
            call require(all(group == [(j + 1, j = 1, size(model%breaks))]), 'left-inclusive cuts')
         end if
         if (m == 1) then
            call require(all(model%bins%iv >= 0.0_dp), 'nonnegative binary IV contributions')
         end if
      end do
   end do
   options = bin_options(separate_special=.true.)
   x(1:4) = [nan, -999.0_dp, inf, -inf]
   y(5) = nan
   call iso_bin(x, y, model, options, [-999.0_dp, -888.0_dp])
   call require(model%status == bin_ok .and. model%n_special == 4, 'four special classes')
   call require(model%membership(5) == 0 .and. size(model%retained) == 399, 'missing-target alignment')
   call apply_bins(model, [-888.0_dp, nan, inf], group, pred)
   call require(group(1) == 0 .and. ieee_is_nan(pred(1)), 'unseen separate special')
   call require(group(2) == 1 .and. group(3) == 3, 'seen nonfinite special classes')
   options%separate_special = .false.
   call iso_bin(x, y, model, options, [-999.0_dp, -888.0_dp])
   call apply_bins(model, [-888.0_dp, nan, inf], group)
   call require(all(group == 1), 'together special prediction')
   ! Invalid inputs never produce a fitted summary.
   call iso_bin(x, y(:399), model)
   call require(model%status == bin_invalid .and. size(model%bins) == 0, 'shape failure')
   call iso_bin(x(:0), y(:0), model)
   call require(model%status == bin_invalid, 'empty input')
   call iso_bin([1.0_dp, 1.0_dp], [0.0_dp, 1.0_dp], model)
   call require(model%status == bin_constant_x, 'constant predictor diagnostic')
   call iso_bin([1.0_dp, 2.0_dp], [0.0_dp, 0.0_dp], model)
   call require(model%status == bin_constant_y, 'constant target diagnostic')
   call iso_bin([nan, inf], [0.0_dp, 1.0_dp], model)
   call require(model%status == bin_no_data, 'no complete observations')
   call iso_bin([1.0_dp, 2.0_dp], [0.0_dp, inf], model)
   call require(model%status == bin_invalid, 'infinite response rejection')
   call fit_bins(x, y, 'bad', model)
   call require(model%status == bin_invalid, 'unknown algorithm rejection')
   options = bin_options(target=1)
   call iso_bin([1.0_dp, 2.0_dp], [0.0_dp, 2.0_dp], model, options)
   call require(model%status == bin_invalid, 'binary target validation')
   options = bin_options(direction=2)
   call iso_bin(x, y, model, options)
   call require(model%status == bin_invalid, 'direction validation')
   options = bin_options(g=-1)
   call pct_bin(x, y, model, options)
   call require(model%status == bin_invalid, 'negative group count')
   options = bin_options(min_pct_obs=-0.1_dp)
   call iso_bin(x, y, model, options)
   call require(model%status == bin_invalid, 'negative minimum observations')
   options = bin_options(p_value=0.0_dp)
   call sts_bin(x, y, model, options)
   call require(model%status == bin_invalid, 'invalid significance threshold')
   options = bin_options(woe_gap=nan)
   call woe_bin(x, y, model, options)
   call require(model%status == bin_invalid, 'NaN option rejection')
   call apply_bins(model, [1.0_dp], group, pred)
   call require(group(1) == 0 .and. ieee_is_nan(pred(1)), 'invalid model prediction')
   ! Exact separation keeps infinite evidence; no undocumented pseudocount.
   call pct_bin([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], model)
   call require(any(.not. ieee_is_finite(model%bins%woe)), 'zero-cell infinite evidence')
   ! Mixed-sign continuous responses explicitly retain log-domain NaN behavior.
   options = bin_options(target=2, g=2, woe_trend=.false.)
   call pct_bin([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [-3.0_dp, -2.0_dp, 1.0_dp, 2.0_dp], model, options)
   call require(any(ieee_is_nan(model%bins%woe)), 'continuous log-domain behavior')
   options%woe_trend = .true.
   call pct_bin([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [-3.0_dp, -2.0_dp, 1.0_dp, 2.0_dp], model, options)
   call require(model%status == bin_ok .and. size(model%bins) == 1, 'NaN evidence continuation')
   options = bin_options(direction=-1, target=2)
   call cum_bin([(real(i, dp), i = 1, 100)], [(real(i, dp), i = 1, 100)], model, options)
   call require(model%status == bin_ok .and. size(model%bins) == 1, 'contradictory cumulative trend correction')
   print *, 'Independent invariant and failure checks:', failures, 'failures'
   if (failures /= 0) error stop 1
contains
   subroutine require(condition, message)
      logical, intent(in) :: condition !! Predicate that must hold for the documented binning contract.
      character(len=*), intent(in) :: message !! Name of the invariant or validation branch under test.
      if (.not. condition) then
         failures = failures + 1
         print *, 'FAIL: ', message
      end if
   end subroutine require
end program test_invariants
