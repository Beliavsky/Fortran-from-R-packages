program test_unit
   use psych
   implicit none
   real(dp) :: table2(2, 2), rho, lower, upper, q, f_value
   real(dp), allocatable :: reversed(:, :), scaled(:, :), scores(:, :), corr(:, :), values(:)
   real(dp), allocatable :: probabilities(:, :), totals(:), difficulty(:)
   integer, allocatable :: counts(:, :)
   real(dp) :: quartiles(3), icc(6), ms(4)
   integer :: failures
   failures = 0

   call check_close(alpha_to_r(0.8_dp, 5), 4.0_dp / 9.0_dp, 1.0e-12_dp, 'alpha_to_r', failures)
   call check_close(r_to_alpha(0.4_dp, 5), 10.0_dp / 13.0_dp, 1.0e-12_dp, 'r_to_alpha', failures)
   call check_close(fisher_r(fisher_z(0.45_dp)), 0.45_dp, 1.0e-12_dp, 'fisher round trip', failures)
   call check_close(logistic_value(0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp), 0.5_dp, 1.0e-12_dp, 'logistic zero', failures)
   call check_close(irt_probability(0.0_dp, 0.0_dp, 1.0_dp), 0.5_dp, 1.0e-12_dp, 'IRT midpoint', failures)
   call check_close(d_to_t(0.5_dp, n=100), 2.5_dp, 1.0e-12_dp, 'd2t', failures)
   call check_close(t_to_d(2.5_dp, n=100), 0.5_dp, 1.0e-12_dp, 't2d', failures)

   table2 = 25.0_dp
   call tetrachoric_table(table2, rho, correction=0.0_dp)
   call check_close(rho, 0.0_dp, 3.0e-4_dp, 'tetrachoric independence', failures)

   call correlation_ci(0.0_dp, 100, 0.95_dp, lower, upper)
   call check_true(lower < 0.0_dp .and. upper > 0.0_dp, 'correlation interval straddles zero', failures)

   block
      real(dp) :: items(3, 2), keys(2)
      items = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 3.0_dp, 2.0_dp], [3, 2])
      keys = [1.0_dp, -1.0_dp]
      call reverse_code_matrix(keys, items, reversed, minima=[1.0_dp, 2.0_dp], maxima=[4.0_dp, 4.0_dp])
      call check_close(reversed(1, 2), 2.0_dp, 1.0e-12_dp, 'reverse code', failures)
      call rescale_matrix(items, 100.0_dp, 15.0_dp, scaled)
      call check_close(sum(scaled(:, 1)) / 3.0_dp, 100.0_dp, 1.0e-10_dp, 'rescale mean', failures)
   end block

   block
      real(dp) :: x(5)
      x = [1.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 3.0_dp]
      q = interpolated_quantile(x, 0.5_dp, 1.0_dp)
      call interpolated_quartiles(x, quartiles, 1.0_dp)
      call check_true(q > 1.5_dp .and. q < 2.5_dp, 'interpolated median', failures)
      call check_close(quartiles(2), q, 1.0e-12_dp, 'interpolated quartile median', failures)
   end block

   block
      real(dp) :: x(4, 2)
      x = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [4, 2])
      call icc_complete_matrix(x, icc, ms)
      call check_close(icc(3), 1.0_dp, 1.0e-12_dp, 'ICC3 perfect agreement', failures)
   end block

   block
      integer :: xcat(8), ycat(8)
      xcat = [1, 1, 1, 1, 2, 2, 2, 2]
      ycat = [1, 1, 2, 2, 1, 1, 2, 2]
      call polychoric_categories(xcat, ycat, rho, correction=0.0_dp)
      call check_close(rho, 0.0_dp, 3.0e-3_dp, 'polychoric independence', failures)
   end block

   block
      real(dp) :: diff(2), disc(2), theta(2)
      diff = [-0.5_dp, 0.5_dp]
      disc = 1.0_dp
      theta = [-1.0_dp, 1.0_dp]
      call irt_expected_scores(theta, diff, disc, probabilities, totals)
      call check_true(totals(2) > totals(1), 'IRT expected score monotonicity', failures)
   end block

   block
      real(dp) :: r(3, 3)
      r = reshape([1.0_dp, 0.5_dp, 0.25_dp, 0.5_dp, 1.0_dp, 0.5_dp, 0.25_dp, 0.5_dp, 1.0_dp], [3, 3])
      call simplex_correlation(3, 0.5_dp, corr)
      call check_matrix(r, corr, 1.0e-12_dp, 'simplex correlation', failures)
   end block

   block
      real(dp) :: x(4, 2), key_matrix(2, 1)
      x = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp], [4, 2])
      key_matrix(:, 1) = [1.0_dp, -1.0_dp]
      call score_key_matrix(key_matrix, x, scores, totals=.false., minima=[1.0_dp, 1.0_dp], maxima=[4.0_dp, 4.0_dp], &
         response_counts=counts)
      call check_close(scores(1, 1), 1.0_dp, 1.0e-12_dp, 'score reverse key', failures)
      call check_true(all(counts(:, 1) == 2), 'score response counts', failures)
   end block

   block
      real(dp) :: item_matrix(4, 2)
      item_matrix = reshape([0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [4, 2])
      call rasch_item_difficulty(item_matrix, difficulty)
      call check_close(difficulty(1), 0.0_dp, 1.0e-12_dp, 'Rasch difficulty 0.5', failures)
      call check_close(difficulty(2), 0.0_dp, 1.0e-12_dp, 'Rasch difficulty second', failures)
   end block

   f_value = 0.0_dp
   values = [f_value]
   if (size(values) /= 1) failures = failures + 1
   if (failures /= 0) error stop 'psych unit tests failed'
   print '(a)', 'All psych unit tests passed'
contains
   subroutine check_close(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual !! Computed scalar under test.
      real(dp), intent(in) :: expected !! Expected reference scalar.
      real(dp), intent(in) :: tolerance !! Maximum accepted absolute difference.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      integer, intent(inout) :: failures !! Accumulated assertion failure count.
      if (abs(actual - expected) > tolerance) then
         print '(a,2es18.8)', trim(label)//' failed: ', actual, expected
         failures = failures + 1
      end if
   end subroutine check_close

   subroutine check_true(condition, label, failures)
      logical, intent(in) :: condition !! Boolean assertion result.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      integer, intent(inout) :: failures !! Accumulated assertion failure count.
      if (.not. condition) then
         print '(a)', trim(label)//' failed'
         failures = failures + 1
      end if
   end subroutine check_true

   subroutine check_matrix(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual(:, :) !! Computed matrix under test.
      real(dp), intent(in) :: expected(:, :) !! Expected reference matrix.
      real(dp), intent(in) :: tolerance !! Maximum accepted elementwise absolute difference.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      integer, intent(inout) :: failures !! Accumulated assertion failure count.
      if (maxval(abs(actual - expected)) > tolerance) then
         print '(a)', trim(label)//' failed'
         failures = failures + 1
      end if
   end subroutine check_matrix
end program test_unit
