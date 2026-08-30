module test_utilities_callbacks
   use learnbayes, only: dp
   implicit none
   private
   public :: gaussian_likelihood
   public :: standard_normal_log_target

contains

   function gaussian_likelihood(observation, parameter1, parameter2, data, params) result(value)
      real(dp), intent(in) :: observation(:) !! Scalar observation stored in element one.
      real(dp), intent(in) :: parameter1 !! Candidate Gaussian mean.
      real(dp), intent(in) :: parameter2 !! Positive candidate Gaussian standard deviation.
      real(dp), intent(in) :: data(:, :) !! Unused numeric data context supplied by the callback API.
      real(dp), intent(in) :: params(:) !! Unused numeric constant context supplied by the callback API.
      real(dp) :: value

      if (size(data) + size(params) < 0) error stop 'unreachable'
      value = exp(-0.5_dp*((observation(1) - parameter1)/parameter2)**2)/parameter2
   end function gaussian_likelihood

   function standard_normal_log_target(theta, data, params) result(value)
      real(dp), intent(in) :: theta(:) !! Parameter vector evaluated under independent standard-normal components.
      real(dp), intent(in) :: data(:, :) !! Unused numeric data context supplied by the callback API.
      real(dp), intent(in) :: params(:) !! Unused numeric constant context supplied by the callback API.
      real(dp) :: value

      if (size(data) + size(params) < 0) error stop 'unreachable'
      value = -0.5_dp*sum(theta*theta)
   end function standard_normal_log_target

end module test_utilities_callbacks

program test_utilities
   use learnbayes
   use test_utilities_callbacks
   implicit none

   type(likelihood_callback) :: likelihood
   type(log_density_callback) :: target
   type(bayes_grid_result) :: grid_result
   real(dp), allocatable :: grouped(:, :)
   integer, allocatable :: player_ids(:)
   integer, allocatable :: seasons(:)
   real(dp), allocatable :: y(:, :)
   real(dp), allocatable :: n_ab(:, :)
   real(dp), allocatable :: age(:, :)
   real(dp) :: data(5, 13)
   integer :: players(5)
   real(dp) :: prior(3, 2)
   real(dp) :: observations(2, 1)
   real(dp) :: xgrid(9)
   real(dp) :: ygrid(9)
   real(dp) :: zgrid(9, 9)
   real(dp) :: p(50)
   real(dp) :: prior_density(50)
   real(dp) :: likelihood_density(50)
   real(dp) :: posterior_density(50)
   integer :: info

   call regroup(reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
      10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp], [5, 2]), 2, grouped)
   call assert_true(size(grouped, 1) == 2 .and. size(grouped, 2) == 2, 'regroup shape')
   call assert_close(sum(grouped), 165.0_dp, 2.0e-14_dp, 'regroup preserves total')

   data = 0.0_dp
   players = [10, 10, 20, 20, 20]
   data(:, 3) = [21.0_dp, 22.0_dp, 24.0_dp, 25.0_dp, 26.0_dp]
   data(:, 5) = [100.0_dp, 110.0_dp, 90.0_dp, 95.0_dp, 105.0_dp]
   data(:, 10) = [20.0_dp, 25.0_dp, 15.0_dp, 18.0_dp, 22.0_dp]
   data(:, 13) = [5.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, 5.0_dp]
   call careertraj_setup(players, data, player_ids, y, n_ab, age, seasons, info)
   call assert_true(info == 0, 'career trajectory setup status')
   call assert_true(all(player_ids == [10, 20]), 'career trajectory IDs')
   call assert_true(all(seasons == [2, 3]), 'career trajectory season counts')

   likelihood%eval => gaussian_likelihood
   observations(:, 1) = [0.2_dp, -0.1_dp]
   prior = 1.0_dp/6.0_dp
   call discrete_bayes_2(likelihood, [-1.0_dp, 0.0_dp, 1.0_dp], [0.5_dp, 1.5_dp], prior, observations, grid_result)
   call assert_close(sum(grid_result%prob), 1.0_dp, 2.0e-14_dp, 'two-parameter Bayes normalization')
   call assert_true(grid_result%predictive > 0.0_dp, 'two-parameter Bayes predictive')

   target%eval => standard_normal_log_target
   call contour_grid(target, [-2.0_dp, 2.0_dp, -2.0_dp, 2.0_dp], 9, xgrid, ygrid, zgrid)
   call assert_close(maxval(zgrid), 0.0_dp, 2.0e-14_dp, 'contour grid normalization')
   call assert_close(xgrid(5), 0.0_dp, 2.0e-14_dp, 'contour x midpoint')
   call assert_close(ygrid(5), 0.0_dp, 2.0e-14_dp, 'contour y midpoint')

   call triplot_data(2.0_dp, 3.0_dp, 4, 6, p, prior_density, likelihood_density, posterior_density)
   call assert_true(all(p > 0.0_dp .and. p < 1.0_dp), 'triplot probability grid')
   call assert_true(all(prior_density >= 0.0_dp), 'triplot prior densities')
   call assert_true(all(likelihood_density >= 0.0_dp), 'triplot likelihood densities')
   call assert_true(all(posterior_density >= 0.0_dp), 'triplot posterior densities')

   print '(a)', 'test_utilities: PASS'

contains

   subroutine assert_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual !! Computed scalar value under test.
      real(dp), intent(in) :: expected !! Reference scalar value expected from the translated algorithm.
      real(dp), intent(in) :: tol !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Short test label reported if the comparison fails.

      if (abs(actual - expected) > tol) then
         write (*, '(a,2es24.15)') trim(label)//' failed: ', actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Boolean condition that must be true for the test to pass.
      character(len=*), intent(in) :: label !! Short test label reported if the condition is false.

      if (.not. condition) then
         write (*, '(a)') trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true

end program test_utilities
