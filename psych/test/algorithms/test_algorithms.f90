program test_algorithms
   use psych
   implicit none
   integer :: failures, info
   type(alpha_result) :: reliability
   type(pca_result) :: pc
   type(cosinor_result) :: cosfit
   type(circadian_result) :: cstats
   real(dp), allocatable :: loadings(:, :), uniqueness(:), theta(:), fit(:)
   real(dp) :: omega_h, omega_total, omega_limit, alpha_value, lambda6
   failures = 0

   block
      real(dp) :: x(5, 3)
      x(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      x(:, 2) = 2.0_dp * x(:, 1)
      x(:, 3) = -x(:, 1)
      call cronbach_alpha_data(x(:, 1:2), reliability, info)
      call check_true(info == 0, 'alpha info', failures)
      call check_close(reliability%standardized_alpha, 1.0_dp, 1.0e-8_dp, 'alpha perfect items', failures)
      call principal_components_data(x, 1, pc)
      call check_true(pc%info == 0, 'PCA info', failures)
      call check_true(pc%values(1) > 2.9_dp, 'PCA first eigenvalue', failures)
   end block

   block
      real(dp) :: r(4, 4)
      integer :: i
      r = 0.4_dp
      do i = 1, 4
         r(i, i) = 1.0_dp
      end do
      allocate(loadings(4, 1), uniqueness(4))
      call principal_axis_factor(r, 1, loadings, uniqueness, info)
      call check_true(info == 0, 'principal-axis factor info', failures)
      call check_true(all(uniqueness >= 0.0_dp), 'principal-axis uniqueness', failures)
      call omega_statistics(r, loadings, omega_h, omega_total, omega_limit, alpha_value, lambda6)
      call check_true(omega_total >= omega_h, 'omega total >= hierarchical', failures)
      call check_true(alpha_value > 0.0_dp, 'omega alpha positive', failures)
      deallocate(loadings, uniqueness)
   end block

   block
      real(dp) :: angle(8), y(8)
      integer :: i
      do i = 1, 8
         angle(i) = real(i - 1, dp) * 3.0_dp
         y(i) = cos(2.0_dp * 3.14159265358979323846_dp * angle(i) / 24.0_dp)
      end do
      call cosinor_fit(angle, y, cosfit, hours=.true., period=24.0_dp)
      call check_true(cosfit%fit > 0.999_dp, 'cosinor exact cosine fit', failures)
      call circadian_stats_vector([23.0_dp, 1.0_dp], cstats, hours=.true., period=24.0_dp)
      call check_true(cstats%mean_phase < 1.0_dp .or. cstats%mean_phase > 23.0_dp, 'circadian wraparound mean', failures)
   end block

   block
      real(dp) :: difficulty(4), responses(3, 4)
      difficulty = [-1.0_dp, -0.3_dp, 0.3_dp, 1.0_dp]
      responses(1, :) = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
      responses(2, :) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
      responses(3, :) = [0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
      call rasch_person_scores(difficulty, responses, theta, fit)
      call check_true(theta(3) > theta(2) .and. theta(2) > theta(1), 'Rasch score ordering', failures)
   end block

   block
      real(dp) :: table2(2, 2), rho
      table2 = reshape([40.0_dp, 10.0_dp, 10.0_dp, 40.0_dp], [2, 2])
      call tetrachoric_table(table2, rho, correction=0.0_dp)
      call check_true(rho > 0.5_dp, 'tetrachoric positive association', failures)
   end block

   if (failures /= 0) error stop 'psych algorithm tests failed'
   print '(a)', 'All psych algorithm tests passed'
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
end program test_algorithms
