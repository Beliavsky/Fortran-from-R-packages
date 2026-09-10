program test_simulation_dgps
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use grf, only : dp, generate_causal_data
   implicit none

   integer, parameter :: n = 96
   integer, parameter :: p = 6
   character(len=10), parameter :: designs(12) = [character(len=10) :: &
      'simple', 'aw1', 'aw2', 'aw3', 'aw3reverse', 'ai1', 'ai2', 'kunzel', &
      'nw1', 'nw2', 'nw3', 'nw4']
   real(dp), allocatable :: x(:,:)
   real(dp), allocatable :: x_copy(:,:)
   real(dp), allocatable :: y(:)
   real(dp), allocatable :: y_copy(:)
   real(dp), allocatable :: tau(:)
   real(dp), allocatable :: tau_copy(:)
   real(dp), allocatable :: baseline(:)
   real(dp), allocatable :: baseline_copy(:)
   real(dp), allocatable :: propensity(:)
   real(dp), allocatable :: propensity_copy(:)
   integer, allocatable :: treatment(:)
   integer, allocatable :: treatment_copy(:)
   integer :: failures
   integer :: info
   integer :: k

   failures = 0
   do k = 1, size(designs)
      call generate_causal_data(n, p, x, y, treatment, tau, baseline, propensity, &
                                info, seed=1200 + k, dgp=trim(designs(k)))
      call check(info == 0, trim(designs(k)) // ' accepted', failures)
      if (info /= 0) cycle
      call check(all(shape(x) == [n,p]), trim(designs(k)) // ' predictor shape', failures)
      call check(all(ieee_is_finite(x)) .and. all(ieee_is_finite(y)), &
                 trim(designs(k)) // ' finite observations', failures)
      call check(all(ieee_is_finite(tau)) .and. all(ieee_is_finite(baseline)), &
                 trim(designs(k)) // ' finite targets', failures)
      call check(all(propensity >= 0.0_dp) .and. all(propensity <= 1.0_dp), &
                 trim(designs(k)) // ' propensity range', failures)
      call check(all(treatment == 0 .or. treatment == 1), &
                 trim(designs(k)) // ' binary treatment', failures)
      call check(abs(sample_sd(baseline) - 1.0_dp) < 1.0e-11_dp, &
                 trim(designs(k)) // ' baseline scaling', failures)
      if (sample_sd(tau) > 1.0e-14_dp) then
         call check(abs(sample_sd(tau) - 0.1_dp) < 1.0e-11_dp, &
                    trim(designs(k)) // ' effect scaling', failures)
      end if
   end do

   call generate_causal_data(n, p, x, y, treatment, tau, baseline, propensity, &
                             info, seed=987, dgp='ai2', sigma_m=1.7_dp, &
                             sigma_tau=0.35_dp, sigma_noise=0.0_dp)
   x_copy = x
   y_copy = y
   treatment_copy = treatment
   tau_copy = tau
   baseline_copy = baseline
   propensity_copy = propensity
   call check(abs(sample_sd(baseline) - 1.7_dp) < 1.0e-11_dp, &
              'custom baseline scaling', failures)
   call check(abs(sample_sd(tau) - 0.35_dp) < 1.0e-11_dp, &
              'custom effect scaling', failures)
   call check(maxval(abs(y - baseline - (real(treatment,dp) - propensity) * tau)) < 1.0e-13_dp, &
              'zero-noise outcome identity', failures)

   call generate_causal_data(n, p, x, y, treatment, tau, baseline, propensity, &
                             info, seed=987, dgp='ai2', sigma_m=1.7_dp, &
                             sigma_tau=0.35_dp, sigma_noise=0.0_dp)
   call check(maxval(abs(x - x_copy)) <= 0.0_dp .and. &
              maxval(abs(y - y_copy)) <= 0.0_dp, &
              'repeatable real-valued simulation', failures)
   call check(all(treatment == treatment_copy) .and. &
              maxval(abs(tau - tau_copy)) <= 0.0_dp, &
              'repeatable treatment and effects', failures)
   call check(maxval(abs(baseline - baseline_copy)) <= 0.0_dp .and. &
              maxval(abs(propensity - propensity_copy)) <= 0.0_dp, &
              'repeatable nuisance components', failures)

   call generate_causal_data(10, 1, x, y, treatment, tau, baseline, propensity, &
                             info, dgp='ai2')
   call check(info == -2, 'minimum predictor count validation', failures)
   call generate_causal_data(10, 3, x, y, treatment, tau, baseline, propensity, &
                             info, dgp='unknown')
   call check(info == -3, 'unknown selector validation', failures)
   call generate_causal_data(10, 3, x, y, treatment, tau, baseline, propensity, &
                             info, sigma_noise=-1.0_dp)
   call check(info == -4, 'negative scale validation', failures)

   if (failures /= 0) error stop 'GRF simulation DGP tests failed'
   print '(a)', 'All GRF causal simulation DGP tests passed.'

contains

   pure real(dp) function sample_sd(values) result(value)
      real(dp), intent(in) :: values(:) !! Values whose sample standard deviation is returned.
      real(dp) :: average

      if (size(values) < 2) then
         value = 0.0_dp
         return
      end if
      average = sum(values) / real(size(values), dp)
      value = sqrt(sum((values - average)**2) / real(size(values) - 1, dp))
   end function sample_sd

   subroutine check(condition, label, failures)
      logical, intent(in) :: condition !! Assertion result.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      integer, intent(inout) :: failures !! Running failure count.

      if (.not. condition) then
         failures = failures + 1
         print '(a,a)', 'FAIL: ', trim(label)
      end if
   end subroutine check

end program test_simulation_dgps
