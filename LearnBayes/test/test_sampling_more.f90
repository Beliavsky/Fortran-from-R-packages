module test_sampling_more_callbacks
   use learnbayes, only: dp
   implicit none
   private
   public :: standard_normal_log_target

contains

   function standard_normal_log_target(theta, data, params) result(value)
      real(dp), intent(in) :: theta(:) !! Parameter vector evaluated under independent standard-normal components.
      real(dp), intent(in) :: data(:, :) !! Unused numeric data context supplied by the callback API.
      real(dp), intent(in) :: params(:) !! Unused numeric constant context supplied by the callback API.
      real(dp) :: value

      if (size(data) + size(params) < 0) error stop 'unreachable'
      value = -0.5_dp*real(size(theta), dp)*log(2.0_dp*acos(-1.0_dp)) - 0.5_dp*sum(theta*theta)
   end function standard_normal_log_target

end module test_sampling_more_callbacks

program test_sampling_more
   use learnbayes
   use test_sampling_more_callbacks
   implicit none

   type(log_density_callback) :: target
   type(mcmc_result) :: chain
   type(rng_state) :: rng
   real(dp) :: proposal_var1(1, 1)
   real(dp) :: proposal_var2(2, 2)
   real(dp), allocatable :: accepted(:, :)
   real(dp) :: simx(2000)
   real(dp) :: simy(2000)
   integer :: info
   integer :: n_accept

   target%eval => standard_normal_log_target
   call rng_seed(rng, 81923_i8)

   call gibbs(rng, target, [0.0_dp, 0.0_dp], 3000, [0.8_dp, 0.8_dp], chain)
   call assert_true(all(chain%accept_by_parameter > 0.2_dp), 'gibbs lower acceptance')
   call assert_true(all(chain%accept_by_parameter < 0.95_dp), 'gibbs upper acceptance')
   call assert_true(abs(sum(chain%par(:, 1))/3000.0_dp) < 0.2_dp, 'gibbs first-coordinate mean')
   call assert_true(abs(sum(chain%par(:, 2))/3000.0_dp) < 0.2_dp, 'gibbs second-coordinate mean')

   proposal_var2 = 0.0_dp
   proposal_var2(1, 1) = 1.0_dp
   proposal_var2(2, 2) = 1.0_dp
   call indepmetrop(rng, target, [0.0_dp, 0.0_dp], proposal_var2, [0.0_dp, 0.0_dp], 1000, chain, info)
   call assert_true(info == 0, 'independence Metropolis status')
   call assert_close(chain%accept_rate, 1.0_dp, 1.0e-14_dp, 'independence Metropolis exact proposal')

   proposal_var1(1, 1) = 1.0_dp
   call rejectsampling(rng, target, [0.0_dp], proposal_var1, 5.0_dp, 0.1_dp, 2000, accepted, n_accept, info)
   call assert_true(info == 0, 'rejection sampler status')
   call assert_true(n_accept > 500 .and. n_accept <= 2000, 'rejection sampler accepted count')
   call assert_true(abs(sum(accepted(:, 1))/real(n_accept, dp)) < 0.2_dp, 'rejection sampler mean')

   call simcontour(rng, target, [-4.0_dp, 4.0_dp, -4.0_dp, 4.0_dp], 2000, simx, simy)
   call assert_true(abs(sum(simx)/2000.0_dp) < 0.2_dp, 'simcontour x mean')
   call assert_true(abs(sum(simy)/2000.0_dp) < 0.2_dp, 'simcontour y mean')
   call assert_true(maxval(abs(simx)) < 4.2_dp, 'simcontour x bounds')
   call assert_true(maxval(abs(simy)) < 4.2_dp, 'simcontour y bounds')

   print '(a)', 'test_sampling_more: PASS'

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

end program test_sampling_more
