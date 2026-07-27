! SPDX-License-Identifier: LGPL-3.0-or-later
program test_hypothesis
   use sharper, only: dp, sr_result, test_result
   use sharper, only: fit_sr, paired_sr_test, unpaired_sr_test
   use sharper, only: sr_equality_test, sr_max_test, sr_conditional_test
   implicit none
   real(dp) :: x(120), y(120), data(120,3), contrast(1,3), covariance(3,3)
   type(sr_result) :: samples(2)
   type(test_result) :: result
   integer :: i

   do i = 1, 120
      x(i) = 0.006_dp+0.018_dp*sin(0.17_dp*real(i,dp))+0.005_dp*cos(0.49_dp*real(i,dp))
      y(i) = 0.004_dp+0.014_dp*sin(0.17_dp*real(i,dp))+0.007_dp*cos(0.31_dp*real(i,dp))
      data(i,1) = x(i)
      data(i,2) = y(i)
      data(i,3) = 0.003_dp+0.01_dp*cos(0.23_dp*real(i,dp))
   end do
   samples(1) = fit_sr(x)
   samples(2) = fit_sr(y)
   result = unpaired_sr_test(samples,alternative='two.sided')
   call assert_probability(result%p_value)
   result = paired_sr_test(x,y,alternative='two.sided')
   call assert_probability(result%p_value)

   contrast(1,:) = [1.0_dp,-1.0_dp,0.0_dp]
   result = sr_equality_test(data,'t','two.sided',contrast)
   call assert_probability(result%p_value)
   result = sr_equality_test(data,'chisq')
   call assert_probability(result%p_value)

   result = sr_max_test([0.5_dp,0.7_dp,0.4_dp],119,method='bonferroni',loglog=.false.)
   call assert_probability(result%p_value)
   result = sr_max_test([0.5_dp,0.7_dp,0.4_dp],119,method='chi-bar-square',loglog=.false.)
   call assert_probability(result%p_value)
   result = sr_max_test([0.5_dp,0.7_dp,0.4_dp],119,method='follman',loglog=.false.)
   call assert_probability(result%p_value)

   covariance = 0.002_dp
   covariance(1,1) = 0.01_dp
   covariance(2,2) = 0.01_dp
   covariance(3,3) = 0.01_dp
   result = sr_conditional_test([0.5_dp,0.7_dp,0.4_dp],119,covariance,alternative='greater')
   call assert_probability(result%p_value)
   if (result%conf_low > result%conf_high) error stop 1

   print '(a)', 'test_hypothesis: PASS'
contains
   subroutine assert_probability(value)
      real(dp), intent(in) :: value
      if (value < 0.0_dp .or. value > 1.0_dp) then
         print '(a,es24.16)', 'invalid probability: ',value
         error stop 1
      end if
   end subroutine assert_probability
end program test_hypothesis
