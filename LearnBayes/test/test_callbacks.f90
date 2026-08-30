module test_callbacks_module
   use learnbayes, only: dp
   implicit none
   private
   public :: bernoulli_likelihood
   public :: normal_log_target
   public :: square_statistic

contains

   function normal_log_target(theta, data, params) result(value)
      real(dp), intent(in) :: theta(:) !! Parameter vector; the first element is evaluated under N(0,1).
      real(dp), intent(in) :: data(:, :) !! Unused numeric data context supplied by the callback API.
      real(dp), intent(in) :: params(:) !! Unused numeric constant context supplied by the callback API.
      real(dp) :: value

      if (size(data) + size(params) < 0) error stop 'unreachable'
      value = -0.5_dp*log(2.0_dp*acos(-1.0_dp)) - 0.5_dp*theta(1)*theta(1)
   end function normal_log_target

   function square_statistic(theta, data, params) result(value)
      real(dp), intent(in) :: theta(:) !! Parameter vector whose first component is squared.
      real(dp), intent(in) :: data(:, :) !! Unused numeric data context supplied by the callback API.
      real(dp), intent(in) :: params(:) !! Unused numeric constant context supplied by the callback API.
      real(dp) :: value

      if (size(data) + size(params) < 0) error stop 'unreachable'
      value = theta(1)*theta(1)
   end function square_statistic

   function bernoulli_likelihood(observation, parameter1, parameter2, data, params) result(value)
      real(dp), intent(in) :: observation(:) !! Scalar Bernoulli observation stored in element one.
      real(dp), intent(in) :: parameter1 !! Candidate Bernoulli success probability.
      real(dp), intent(in) :: parameter2 !! Unused second parameter supplied by the generic discrete-Bayes API.
      real(dp), intent(in) :: data(:, :) !! Unused numeric data context supplied by the callback API.
      real(dp), intent(in) :: params(:) !! Unused numeric constant context supplied by the callback API.
      real(dp) :: value

      if (parameter2 + real(size(data) + size(params), dp) < -huge(1.0_dp)) error stop 'unreachable'
      if (observation(1) > 0.5_dp) then
         value = parameter1
      else
         value = 1.0_dp - parameter1
      end if
   end function bernoulli_likelihood

end module test_callbacks_module

program test_callbacks
   use learnbayes
   use test_callbacks_module
   implicit none

   type(log_density_callback) :: target
   type(log_density_callback) :: stat
   type(likelihood_callback) :: like
   type(laplace_result) :: lap
   type(mcmc_result) :: chain
   type(importance_result) :: imp
   type(bayes_discrete_result) :: db
   type(rng_state) :: rng
   real(dp) :: proposal_var(1, 1)
   real(dp) :: theta_sir(1000, 1)
   real(dp) :: observations(3, 1)
   real(dp) :: parameter(3)
   real(dp) :: prior(3)
   integer :: info

   call rng_seed(rng, 1234567_i8)
   target%eval => normal_log_target
   stat%eval => square_statistic
   call laplace(target, [0.8_dp], lap)
   call assert_true(lap%converged, 'laplace convergence')
   call assert_close(lap%mode(1), 0.0_dp, 2.0e-5_dp, 'laplace mode')
   call assert_close(lap%var(1, 1), 1.0_dp, 2.0e-5_dp, 'laplace variance')
   call assert_close(lap%log_integral, 0.0_dp, 3.0e-5_dp, 'laplace normalized integral')

   proposal_var(1, 1) = 1.0_dp
   call rwmetrop(rng, target, proposal_var, 1.0_dp, [0.0_dp], 2000, chain, info)
   call assert_true(info == 0, 'rwmetrop status')
   call assert_true(chain%accept_rate > 0.2_dp .and. chain%accept_rate < 0.95_dp, 'rwmetrop acceptance')
   call assert_true(abs(sum(chain%par(:, 1))/real(size(chain%par, 1), dp)) < 0.25_dp, 'rwmetrop mean')

   call impsampling(rng, target, stat, [0.0_dp], proposal_var, 5.0_dp, 5000, imp, info)
   call assert_true(info == 0, 'importance status')
   call assert_true(abs(imp%estimate - 1.0_dp) < 0.12_dp, 'importance E(theta^2)')

   call sir(rng, target, [0.0_dp], proposal_var, 5.0_dp, 1000, theta_sir, info)
   call assert_true(info == 0, 'sir status')
   call assert_true(abs(sum(theta_sir(:, 1))/1000.0_dp) < 0.25_dp, 'sir mean')

   like%eval => bernoulli_likelihood
   observations(:, 1) = [1.0_dp, 1.0_dp, 0.0_dp]
   parameter = [0.25_dp, 0.5_dp, 0.75_dp]
   prior = 1.0_dp/3.0_dp
   call discrete_bayes(like, parameter, prior, observations, db)
   call assert_close(sum(db%prob), 1.0_dp, 2.0e-14_dp, 'discrete bayes normalization')
   call assert_true(db%prob(3) > db%prob(1), 'discrete bayes ordering')

   print '(a)', 'test_callbacks: PASS'

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

end program test_callbacks
