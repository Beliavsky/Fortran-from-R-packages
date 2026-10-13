! SPDX-License-Identifier: GPL-2.0-or-later
program test_outliers
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use outliers_api
   implicit none
   integer :: failures, checks
   failures = 0
   checks = 0
   call test_distributions()
   call test_wrappers()
   call test_scores()
   call test_helpers()
   print '(a,i0,a,i0)', 'outliers: ', checks, ' checks, failures: ', failures
   if (failures /= 0) error stop 1
contains
   subroutine close_value(got, wanted, tolerance, label)
      real(dp), intent(in) :: got !! Computed numerical result.
      real(dp), intent(in) :: wanted !! Actual-R fixture value, including NaN or infinity.
      real(dp), intent(in) :: tolerance !! Absolute comparison tolerance for this algorithm.
      character(len = *), intent(in) :: label !! Failure label identifying the numerical operation.
      logical :: ok
      checks = checks + 1
      if (ieee_is_nan(wanted)) then
         ok = ieee_is_nan(got)
      else if (.not. ieee_is_finite(wanted)) then
         ok = got == wanted
      else
         ok = ieee_is_finite(got)
         if (ok) ok = abs(got - wanted) <= tolerance
      end if
      if (.not. ok) then
         failures = failures + 1
         if (failures <= 100) print '(a,2es25.16)', trim(label) // ': ', got, wanted
      end if
   end subroutine close_value

   subroutine require(ok, label)
      logical, intent(in) :: ok !! Logical property expected to hold.
      character(len = *), intent(in) :: label !! Failure description.
      checks = checks + 1
      if (.not. ok) then
         failures = failures + 1
         if (failures <= 100) print '(a)', label
      end if
   end subroutine require

   subroutine test_distributions()
      integer :: unit, ios, code, n, kind, k, rows
      real(dp) :: a, b, wanted, got, tolerance, p, lower, boundary_tolerance
      character(len = 100) :: label
      open(newunit = unit, file = 'test/distributions.dat', status = 'old', action = 'read')
      rows = 0
      do
         read(unit, *, iostat = ios) code, n, kind, k, a, b, wanted
         if (ios < 0) exit
         if (ios /= 0) error stop 'Invalid distribution fixture'
         rows = rows + 1
         tolerance = 2.0e-11_dp
         select case (code)
         case (1)
            got = qdixon(a, n, kind)
            tolerance = 2.0e-10_dp
         case (2)
            got = pdixon(a, n, kind)
            tolerance = 5.0e-8_dp
         case (3)
            got = qgrubbs(a, n, kind)
         case (4)
            got = pgrubbs(a, n, kind)
            if (kind == 11) tolerance = 8.0e-5_dp
         case (5)
            got = qcochran(a, b, k)
         case (6)
            got = pcochran(a, b, k)
         end select
         write(label, '(a,4i6,2es15.5)') 'Distribution ', code, n, kind, k, a, b
         ! At the source lower clipping discontinuity, R itself can yield either
         ! zero or a bracket error depending on one-ulp qt rounding. Both are
         ! explicitly checked; ordinary interior/domain fixtures remain strict.
         if (code == 4 .and. kind == 11) then
            lower = qgrubbs(2.0e-16_dp, n, 11)
            boundary_tolerance = 128.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(a))
            if (abs(a - lower) <= boundary_tolerance) then
               call require(ieee_is_nan(got) .or. got == 0.0_dp, 'Grubbs lower discontinuity result')
               call require(ieee_is_nan(wanted) .or. wanted == 0.0_dp, 'R lower discontinuity fixture')
               cycle
            end if
         end if
         call close_value(got, wanted, tolerance, label)
      end do
      close(unit)
      call require(rows > 10000, 'Exhaustive distribution fixtures loaded')
      do n = 4, 30
         p = 0.925_dp
         call close_value(pgrubbs(qgrubbs(p, n), n), p, 2.0e-11_dp, 'Grubbs 10 inverse')
         call close_value(pgrubbs(qgrubbs(p, n, 11), n, 11), p, 2.0e-11_dp, 'Grubbs 11 inverse')
         call close_value(pcochran(qcochran(p, real(n, dp), 5), real(n, dp), 5), p, &
            2.0e-11_dp, 'Cochran inverse')
      end do
      call close_value(qdixon(0.05_dp, 3), 0.941_dp, 1.0e-14_dp, 'Published Dixon table node')
      call close_value(qdixon(0.7_dp, 7, 21), 1.0_dp, 1.0e-14_dp, 'Preserve source 350 anomaly and clamp')
      call close_value(qdixon(0.4_dp, 6, 22), 1.0_dp, 1.0e-14_dp, 'Preserve source 830 anomaly and clamp')
      call close_value(pdixon(0.0_dp, 29, 12), 0.725_dp, 1.0e-14_dp, 'Tied endpoint regression')
      call require(ieee_is_nan(qdixon(0.05_dp, 31)), 'Dixon upper size restriction')
      call require(ieee_is_nan(qdixon(0.05_dp, 5, 22)), 'Dixon type-specific minimum size')
      call require(ieee_is_nan(qgrubbs(0.05_dp, 31, 20)), 'Grubbs table upper size restriction')
      call require(ieee_is_nan(qgrubbs(0.05_dp, 3, 20)), 'Grubbs table missing size 3')
      call require(ieee_is_nan(qcochran(-0.1_dp, 10.0_dp, 5)), 'Cochran invalid probability')
   end subroutine test_distributions

   subroutine test_wrappers()
      integer :: unit, ios, code, n, kind, opp, two
      real(dp) :: variance, wanted(3), tolerance
      real(dp), allocatable :: x(:)
      type(outlier_test_result) :: res
      character(len = 100) :: label
      open(newunit = unit, file = 'test/tests.dat', status = 'old', action = 'read')
      do
         read(unit, *, iostat = ios) code, n, kind, opp, two, variance, wanted
         if (ios < 0) exit
         if (ios /= 0) error stop 'Invalid test fixture'
         allocate(x(n))
         read(unit, *) x
         select case (code)
         case (1)
            res = grubbs_test(x, kind, opp == 1, two == 1)
         case (2)
            res = dixon_test(x, kind, opp == 1, two == 1)
         case (3)
            if (variance == 0.0_dp) then
               res = chisq_out_test(x, opposite = opp == 1)
            else
               res = chisq_out_test(x, variance, opp == 1)
            end if
         end select
         write(label, '(a,5i4)') 'Test ', code, n, kind, opp, two
         ! These documented invalid-domain choices intentionally replace R's
         ! finite n=2 G and spurious constant-sample p=0 with an invalid result.
         if (code == 1) then
            if (n < 3 .or. all(x - x(1) == 0.0_dp)) then
               wanted = ieee_value(0.0_dp, ieee_quiet_nan)
               call require(res%status == 1, 'Invalid Grubbs size/constant status')
            end if
         end if
         call close_value(res%statistic, wanted(1), 2.0e-11_dp, label // ' statistic')
         if (code == 1) call close_value(res%u, wanted(2), 2.0e-11_dp, label // ' U')
         tolerance = 5.0e-8_dp
         if (code == 1 .and. kind == 11) tolerance = 1.6e-4_dp
         call close_value(res%p_value, wanted(3), tolerance, label // ' p-value')
         deallocate(x)
      end do
      close(unit)
   end subroutine test_wrappers

   subroutine test_scores()
      integer :: unit, ios, n, kind, i
      character(len = 5), parameter :: names(5) = [character(len = 5) :: 'z', 't', 'chisq', 'iqr', 'mad']
      real(dp) :: p, lim
      real(dp), allocatable :: x(:), wanted(:), got(:)
      character(len = 100) :: label
      open(newunit = unit, file = 'test/scores.dat', status = 'old', action = 'read')
      do
         read(unit, *, iostat = ios) n, kind, p, lim
         if (ios < 0) exit
         if (ios /= 0) error stop 'Invalid score fixture'
         allocate(x(n), wanted(n))
         read(unit, *) x
         read(unit, *) wanted
         got = scores(x, type = trim(names(kind)), prob = p, lim = lim)
         write(label, '(a,2i4,2es12.3)') 'Scores ', n, kind, p, lim
         do i = 1, n
            call close_value(got(i), wanted(i), 2.0e-11_dp, label)
         end do
         deallocate(x, wanted, got)
      end do
      close(unit)
   end subroutine test_scores

   subroutine test_helpers()
      real(dp) :: x(6), matrix(6, 2), nan, infinity, got(6, 2), values(2), table(5), probabilities(5)
      real(dp), allocatable :: removed(:)
      logical :: mask(6), mask2(6, 2)
      type(outlier_column) :: columns(2)
      type(outlier_test_result) :: res, other
      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      infinity = ieee_value(0.0_dp, ieee_positive_inf)
      x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 20.0_dp, 20.0_dp]
      call close_value(outlier(x), 20.0_dp, 0.0_dp, 'Outlier selects upper tie')
      call close_value(outlier(x, .true.), 1.0_dp, 0.0_dp, 'Opposite outlier')
      mask = outlier_mask(x)
      call require(all(mask .eqv. [.false., .false., .false., .false., .true., .true.]), 'Mask selects all ties')
      removed = rm_outlier(x)
      call require(size(removed) == 4, 'All tied outliers removed')
      call require(all(removed == x(:4)), 'Removal preserves input order')
      removed = rm_outlier(x, fill = .true.)
      call require(all(removed(5:) == 2.5_dp), 'Mean replacement excludes every tied outlier')
      removed = rm_outlier(x, fill = .true., median = .true., opposite = .true.)
      call close_value(removed(1), 4.0_dp, 0.0_dp, 'Opposite median replacement')
      matrix(:, 1) = x
      matrix(:, 2) = -x
      values = outlier(matrix)
      call require(all(values == [20.0_dp, -20.0_dp]), 'Matrix outlier by column')
      mask2 = outlier_mask(matrix)
      call require(all(mask2(:, 1) .eqv. mask2(:, 2)), 'Matrix tied masks')
      columns = rm_outlier(matrix)
      call require(size(columns(1)%values) == 4 .and. size(columns(2)%values) == 4, 'Ragged column removal')
      got = scores(matrix)
      call require(all(abs(got(:, 1) + got(:, 2)) < 1.0e-14_dp), 'Matrix score sign symmetry')
      x(1) = nan
      call close_value(outlier(x), 20.0_dp, 0.0_dp, 'Missing values omitted for selection')
      mask = outlier_mask(x)
      call require(.not. mask(1), 'Missing logical selection is false')
      removed = rm_outlier(x)
      call require(ieee_is_nan(removed(1)), 'Removal preserves missing observations')
      removed = rm_outlier(x, fill = .true.)
      call require(all(ieee_is_nan(removed(5:))), 'Fill propagates missing remaining values')
      call require(ieee_is_nan(outlier([nan, nan])), 'All-missing endpoint')
      call require(ieee_is_nan(outlier([1.0_dp, infinity])), 'Infinite endpoint is explicitly unsupported')
      removed = rm_outlier([2.0_dp, 2.0_dp])
      call require(size(removed) == 0, 'Constant sample removes all ties')
      call require(size(outlier_mask([real(dp) ::])) == 0, 'Empty sample mask')
      probabilities = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
      table = probabilities ** 3 - 2.0_dp * probabilities + 3.0_dp
      call close_value(qtable(0.4_dp, probabilities, table), 2.264_dp, 2.0e-14_dp, 'Independent cubic interpolation')
      call close_value(qtable(-1.0_dp, probabilities, table), 4.9375_dp, 2.0e-14_dp, 'Linear lower extrapolation')
      call close_value(qtable(0.4_dp, probabilities(5:1:-1), table(5:1:-1)), &
         2.264_dp, 2.0e-14_dp, 'Unsorted table is sorted in pairs')
      call require(all(abs(qtable([0.4_dp, 0.5_dp], probabilities, table) - &
         [2.264_dp, 2.125_dp]) < 2.0e-14_dp), 'Vector interpolation')
      call require(ieee_is_nan(qtable(0.5_dp, [0.0_dp, 1.0_dp], [2.0_dp, 3.0_dp])), &
         'Interior cubic needs four coordinates')
      res = cochran_test([1.0_dp, 2.0_dp, 9.0_dp], [5.0_dp, 5.0_dp, 5.0_dp])
      call close_value(res%statistic, 0.75_dp, 0.0_dp, 'Cochran maximum variance fraction')
      call close_value(res%p_value, 1.0_dp - pcochran(0.75_dp, 5.0_dp, 3), 1.0e-14_dp, 'Cochran upper p-value')
      call require(res%group == 3 .and. res%status == 0, 'Cochran group selection')
      res = cochran_test([1.0_dp, 1.0_dp, 9.0_dp], [5.0_dp, 5.0_dp, 5.0_dp], inlying = .true.)
      call require(all(res%groups == [1, 2]), 'Cochran tied minimum groups')
      call close_value(res%statistic, 1.0_dp / 11.0_dp, 0.0_dp, 'Cochran inlying fraction')
      res = cochran_test(matrix)
      other = cochran_groups(reshape(matrix, [12]), [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2])
      call close_value(res%statistic, other%statistic, 0.0_dp, 'Cochran matrix and grouped agreement')
      call close_value(res%df, 6.0_dp, 0.0_dp, 'Cochran group sample size')
      res = cochran_groups([1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 5.0_dp], [1, 1, 2, 2, 2])
      call close_value(res%df, 2.5_dp, 0.0_dp, 'Unequal group counts averaged')
      res = cochran_groups([1.0_dp, 2.0_dp], [1, 3])
      call require(res%status /= 0, 'Cochran invalid groups rejected')
   end subroutine test_helpers
end program test_outliers
