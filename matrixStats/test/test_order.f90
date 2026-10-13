program test_matrixstats_order
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, &
      ieee_is_nan, ieee_is_finite
   use r_kinds, only: dp
   use matrixstats_order
   implicit none
   real(dp) :: nan, inf, samples(9, 4), probs(10), xs(6, 6), ws(6, 8), ranks(9), expected
   real(dp), allocatable :: empty(:), tie_x(:), tie_w(:)
   integer :: s, t, p, i, na, method, reference, status, comparisons, failures
   character(len=12), parameter :: rank_methods(6) = [character(len=12) :: &
      'max', 'min', 'average', 'first', 'last', 'dense']
   character(len=8), parameter :: median_methods(4) = [character(len=8) :: 'weighted', 'min', 'max', 'mean']
   logical :: remove

   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   comparisons = 0
   failures = 0
   open(newunit=reference, file='test/order_reference.txt', status='old', action='read', iostat=status)
   if (status /= 0) then
      open(newunit=reference, file='order_reference.txt', status='old', action='read', iostat=status)
   end if
   if (status /= 0) error stop 'Cannot read committed R ordering reference fixture'

   samples(:, 1) = [8.0_dp, -2.0_dp, 4.0_dp, 4.0_dp, 9.0_dp, 1.0_dp, 7.0_dp, 0.0_dp, 3.0_dp]
   samples(:, 2) = [nan, 8.0_dp, 1.0_dp, nan, -1.0_dp, 6.0_dp, 6.0_dp, 2.0_dp, 3.0_dp]
   samples(:, 3) = [inf, 3.0_dp, -inf, 2.0_dp, 8.0_dp, 8.0_dp, 0.0_dp, 9.0_dp, 1.0_dp]
   samples(:, 4) = nan
   probs = [0.0_dp, 0.01_dp, 0.1_dp, 0.25_dp, 1.0_dp / 3.0_dp, 0.5_dp, 0.7_dp, 0.9_dp, 0.99_dp, 1.0_dp]
   do s = 1, size(samples, 2)
      do t = 1, 9
         do p = 1, size(probs)
            call check_r(key('q', s, t, p), ms_quantile(samples(:, s), probs(p), na_rm=.true., type=t))
         end do
         call check_r(key('iqr', s, t), ms_iqr(samples(:, s), na_rm=.true., type=t))
      end do
      do na = 0, 1
         remove = na == 1
         call check_r(key('median', s, na), ms_median(samples(:, s), na_rm=remove))
         call check_r(key('mad', s, na), ms_mad(samples(:, s), na_rm=remove))
         call check_r(key('mad_center', s, na), ms_mad(samples(:, s), center=2.0_dp, constant=2.0_dp, na_rm=remove))
      end do
      do method = 1, size(rank_methods)
         ranks = ms_ranks(samples(:, s), ties_method=trim(rank_methods(method)))
         do i = 1, size(samples, 1)
            call check_r(key('rank_' // trim(rank_methods(method)), s, i), ranks(i))
         end do
      end do
      if (.not. any(ieee_is_nan(samples(:, s)))) then
         do i = 1, size(samples, 1)
            call check_r(key('order', s, i), ms_order_stat(samples(:, s), i))
         end do
      end if
   end do

   xs(:, 1) = [-4.0_dp, 7.0_dp, 2.0_dp, 9.0_dp, 2.0_dp, 1.0_dp]
   xs(:, 2) = [nan, 7.0_dp, 2.0_dp, 9.0_dp, 2.0_dp, 1.0_dp]
   xs(:, 3) = [nan, 7.0_dp, 2.0_dp, nan, 2.0_dp, 1.0_dp]
   xs(:, 4) = [inf, 7.0_dp, 2.0_dp, 9.0_dp, 2.0_dp, 1.0_dp]
   xs(:, 5) = [-inf, 7.0_dp, 2.0_dp, inf, 2.0_dp, 1.0_dp]
   xs(:, 6) = nan
   ws(:, 1) = 1.0_dp
   ws(:, 2) = [0.0_dp, 2.0_dp, 4.0_dp, 1.0_dp, 0.0_dp, 3.0_dp]
   ws(:, 3) = [-2.0_dp, 2.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 3.0_dp]
   ws(:, 4) = [nan, 2.0_dp, 4.0_dp, 1.0_dp, 0.0_dp, 3.0_dp]
   ws(:, 5) = [inf, 2.0_dp, 4.0_dp, inf, 0.0_dp, 3.0_dp]
   ws(:, 6) = 0.1_dp
   ws(:, 7) = 0.0_dp
   ws(:, 8) = [inf, 2.0_dp, -inf, 1.0_dp, 0.0_dp, 3.0_dp]
   do s = 1, size(xs, 2)
      do t = 1, size(ws, 2)
         do na = 0, 1
            remove = na == 1
            call check_r(key('wmean', s, t, na), ms_weighted_mean(xs(:, s), ws(:, t), na_rm=remove))
            call check_r(key('wrefine', s, t, na), ms_weighted_mean(xs(:, s), ws(:, t), na_rm=remove, refine=.true.))
            call check_r(key('wmedian', s, t, na), ms_weighted_median(xs(:, s), ws(:, t), na_rm=remove))
            do method = 1, size(median_methods)
               call check_r(key('wmedian_' // trim(median_methods(method)), s, t, na), &
                  ms_weighted_median(xs(:, s), ws(:, t), na_rm=remove, ties=trim(median_methods(method))))
            end do
            call check_r(key('wvar', s, t, na), ms_weighted_var(xs(:, s), ws(:, t), na_rm=remove))
            call check_r(key('wsd', s, t, na), ms_weighted_sd(xs(:, s), ws(:, t), na_rm=remove))
            call check_r(key('wmad', s, t, na), ms_weighted_mad(xs(:, s), ws(:, t), na_rm=remove))
            call check_r(key('wvar_center', s, t, na), &
               ms_weighted_var(xs(:, s), ws(:, t), na_rm=remove, center=2.0_dp))
            call check_r(key('wmad_center', s, t, na), &
               ms_weighted_mad(xs(:, s), ws(:, t), na_rm=remove, center=2.0_dp, constant=2.0_dp))
         end do
      end do
   end do
   do s = 1, size(xs, 2)
      do na = 0, 1
         remove = na == 1
         call check_r(key('unit_mean', s, na), ms_weighted_mean(xs(:, s), na_rm=remove))
         call check_r(key('unit_median', s, na), ms_weighted_median(xs(:, s), na_rm=remove))
         call check_r(key('unit_var', s, na), ms_weighted_var(xs(:, s), na_rm=remove))
         call check_r(key('unit_mad', s, na), ms_weighted_mad(xs(:, s), na_rm=remove))
      end do
   end do
   do s = 1, 80
      tie_x = [(real(mod(7 * i + 3 * s, 11) - 5, dp), i = 1, s)]
      tie_w = [(real(mod(3 * i + 2 * s, 7), dp), i = 1, s)]
      call check_r(key('tie_median', s), ms_weighted_median(tie_x, tie_w))
      call check_r(key('tie_mad', s), ms_weighted_mad(tie_x, tie_w))
      do method = 1, size(median_methods)
         call check_r(key('tie_' // trim(median_methods(method)), s), &
            ms_weighted_median(tie_x, tie_w, ties=trim(median_methods(method))))
      end do
   end do
   call check_r('unequal_tie_mad', &
      ms_weighted_mad([4.0_dp, 1.0_dp, 7.0_dp, 1.0_dp, 9.0_dp, 2.0_dp], &
                      [1.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 4.0_dp, 2.0_dp]))
   read(reference, *, iostat=status) expected
   if (status >= 0) error stop 'Unconsumed R ordering references'
   close(reference)

   ! Independent small examples and invalid-input/special-value contracts.
   allocate(empty(0))
   call check('median even', ms_median([1.0_dp, 4.0_dp, 2.0_dp, 3.0_dp]), 2.5_dp)
   call check('median odd', ms_median([9.0_dp, 3.0_dp, 1.0_dp]), 3.0_dp)
   call check('median empty', ms_median(empty), nan)
   call check('median missing', ms_median([1.0_dp, nan]), nan)
   call check('median remove', ms_median([1.0_dp, nan], .true.), 1.0_dp)
   call check('median equal infinite', ms_median([inf, inf]), inf)
   call check('quantile exact endpoint', ms_quantile([-inf, 1.0_dp], 0.0_dp), -inf)
   call check('quantile equal infinite', ms_quantile([inf, inf], 0.3_dp), inf)
   call check('quantile empty', ms_quantile(empty, 0.5_dp), nan)
   call check('quantile invalid p', ms_quantile([1.0_dp], -0.1_dp), nan)
   call check('quantile missing p', ms_quantile([1.0_dp], nan), nan)
   call check('quantile invalid type', ms_quantile([1.0_dp], 0.5_dp, type=10), nan)
   call check('quantile missing observation', ms_quantile([1.0_dp, nan], 0.5_dp), nan)
   call check('iqr single infinite', ms_iqr([inf]), 0.0_dp)
   call check('iqr independent', ms_iqr([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]), 2.0_dp)
   call check('mad independent', ms_mad([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]), 1.4826_dp)
   call check('mad singleton infinite', ms_mad([inf]), 0.0_dp)
   call check('mad singleton center', ms_mad([2.0_dp], center=0.0_dp, constant=1.0_dp), 2.0_dp)
   call check('mad infinity residual', ms_mad([inf, inf, 1.0_dp], na_rm=.true.), nan)
   call check('order low invalid', ms_order_stat([1.0_dp], 0), nan)
   call check('order high invalid', ms_order_stat([1.0_dp], 2), nan)
   call check('order missing', ms_order_stat([1.0_dp, nan], 1), nan)
   call check('weights mean negative', ms_weighted_mean([1.0_dp, 2.0_dp], [-1.0_dp, 2.0_dp]), 3.0_dp)
   call check('weights mean zero omission', ms_weighted_mean([nan, 2.0_dp], [0.0_dp, 1.0_dp]), 2.0_dp)
   call check('weights mean all zero', ms_weighted_mean([1.0_dp, 2.0_dp], [0.0_dp, 0.0_dp]), nan)
   call check('weights mean zero denominator', ms_weighted_mean([1.0_dp, 2.0_dp], [-1.0_dp, 1.0_dp]), inf)
   call check('weights mean mismatch', ms_weighted_mean([1.0_dp], [1.0_dp, 1.0_dp]), nan)
   call check('weights median mismatch', ms_weighted_median([1.0_dp], [1.0_dp, 1.0_dp]), nan)
   call check('weights variance mismatch', ms_weighted_var([1.0_dp], [1.0_dp, 1.0_dp]), nan)
   call check('weights median interpolate', ms_weighted_median([0.0_dp, 10.0_dp], [1.0_dp, 3.0_dp]), 7.5_dp)
   call check('weights median classical', &
      ms_weighted_median([0.0_dp, 10.0_dp], [1.0_dp, 3.0_dp], interpolate=.false.), 10.0_dp)
   call check('weights median lower', ms_weighted_median([0.0_dp, 10.0_dp], ties='min'), 0.0_dp)
   call check('weights median upper', ms_weighted_median([0.0_dp, 10.0_dp], ties='max'), 10.0_dp)
   call check('weights median ties override', &
      ms_weighted_median([0.0_dp, 10.0_dp], [1.0_dp, 3.0_dp], interpolate=.true., ties='max'), 7.5_dp)
   call check('weights median invalid ties', ms_weighted_median([1.0_dp], ties='random'), nan)
   call check('weights variance independent', ms_weighted_var([1.0_dp, 3.0_dp], [2.0_dp, 2.0_dp]), 4.0_dp / 3.0_dp)
   call check('weights variance tiny total', ms_weighted_var([0.0_dp, 2.0_dp], [0.25_dp, 0.25_dp]), -1.0_dp)
   call check('weights sd negative variance', ms_weighted_sd([0.0_dp, 2.0_dp], [0.25_dp, 0.25_dp]), nan)
   call check('weights variance unit total', ms_weighted_var([0.0_dp, 2.0_dp], [0.5_dp, 0.5_dp]), inf)
   call check('weights variance constant unit', ms_weighted_var([2.0_dp, 2.0_dp], [0.5_dp, 0.5_dp]), nan)
   call check('weights mad single center', ms_weighted_mad([2.0_dp], [1.0_dp], center=0.0_dp), 0.0_dp)
   call check('weights mad absent center', ms_weighted_mad([2.0_dp], center=0.0_dp, constant=1.0_dp), 2.0_dp)
   call check('weights mad absent infinite', ms_weighted_mad([inf]), nan)
   call check('weights mad robust even median', ms_weighted_mad([1.0e308_dp, 1.0e308_dp]), 0.0_dp)
   call check('weights mad robust residual', &
      ms_weighted_mad([1.0e308_dp, 1.0e308_dp], center=0.0_dp, constant=1.0_dp), 1.0e308_dp)
   ! Deliberate repair of upstream's odd-count infinite-weight partial-sort bug.
   call check('weights infinite odd true median', ms_weighted_median([2.0_dp, 0.0_dp, 2.0_dp], [inf, inf, inf]), 2.0_dp)
   call check('weights infinite odd distinct', ms_weighted_median([1.0_dp, 0.0_dp, 2.0_dp], [inf, inf, inf]), 1.0_dp)
   call check_rank_extensions()
   if (failures > 0) then
      print '(i0, a, i0)', failures, ' ordering/weight checks failed out of ', comparisons
      error stop 1
   end if
   print '(a, i0, a)', 'matrixStats ordering/weight tests: ', comparisons, ' checks passed'
contains
   function key(prefix, i, j, k) result(label)
      character(len=*), intent(in) :: prefix !! Reference family identifier used by the R fixture generator.
      integer, intent(in) :: i !! First deterministic case index.
      integer, intent(in), optional :: j !! Optional second deterministic case index.
      integer, intent(in), optional :: k !! Optional third deterministic case index.
      character(len=96) :: label
      if (present(k)) then
         write(label, '(a,3("_",i0))') prefix, i, j, k
      else if (present(j)) then
         write(label, '(a,2("_",i0))') prefix, i, j
      else
         write(label, '(a,"_",i0)') prefix, i
      end if
   end function key

   subroutine check_r(label, actual)
      character(len=*), intent(in) :: label !! Expected ordered label in the committed R fixture.
      real(dp), intent(in) :: actual !! Fortran result to compare with the next fixture value.
      character(len=96) :: recorded
      integer :: io
      real(dp) :: oracle
      read(reference, *, iostat=io) recorded, oracle
      if (io /= 0) error stop 'Missing or invalid R ordering reference record'
      if (trim(recorded) /= trim(label)) then
         print *, 'Reference order mismatch: ', trim(recorded), ' versus ', trim(label)
         error stop 2
      end if
      call check(label, actual, oracle)
   end subroutine check_r

   subroutine check(label, actual, oracle)
      character(len=*), intent(in) :: label !! Short diagnostic identifying the asserted behavior.
      real(dp), intent(in) :: actual !! Computed numerical result, including permitted NaNs and infinities.
      real(dp), intent(in) :: oracle !! Expected result, with all NaNs treated as equivalent missing values.
      logical :: matches
      comparisons = comparisons + 1
      if (ieee_is_nan(oracle)) then
         matches = ieee_is_nan(actual)
      else if (.not. ieee_is_finite(oracle)) then
         matches = actual == oracle
      else if (.not. ieee_is_finite(actual)) then
         matches = .false.
      else
         matches = abs(actual / max(1.0_dp, abs(oracle)) - oracle / max(1.0_dp, abs(oracle))) <= 2.0e-12_dp
      end if
      if (.not. matches) then
         failures = failures + 1
         print '(a, a, a, es24.16, a, es24.16)', 'FAIL ', trim(label), ': actual=', actual, ' expected=', oracle
      end if
   end subroutine check

   subroutine check_rank_extensions()
      real(dp) :: result(4)
      integer :: j
      real(dp), parameter :: first_expected(4) = [1.0_dp, 4.0_dp, 2.0_dp, 3.0_dp]
      result = ms_ranks([nan, 2.0_dp, nan, 1.0_dp], na_last='first')
      do j = 1, 4
         call check(key('rank missing first', j), result(j), first_expected(j))
      end do
      result = ms_ranks([nan, 2.0_dp, nan, 1.0_dp], na_last='last')
      call check('rank missing last 1', result(1), 3.0_dp)
      call check('rank missing last 2', result(2), 2.0_dp)
      call check('rank missing last 3', result(3), 4.0_dp)
      call check('rank missing last 4', result(4), 1.0_dp)
      result = ms_ranks([1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], ties_method='random')
      do j = 1, 4
         call check(key('rank random unsupported', j), result(j), nan)
      end do
      result = ms_ranks([1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], na_last='invalid')
      do j = 1, 4
         call check(key('rank invalid placement', j), result(j), nan)
      end do
   end subroutine check_rank_extensions
end program test_matrixstats_order
