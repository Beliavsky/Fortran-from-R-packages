! SPDX-License-Identifier: GPL-3.0-only
program test_semiparametric
   use ufrisk
   use test_support
   implicit none
   integer, parameter :: n = 125
   real(dp) :: prices(n),return_value,scale_value
   integer :: i
   type(ufrisk_options) :: options
   type(ufrisk_result) :: short_result,long_result

   prices(1) = 100.0_dp
   do i = 2,n
      scale_value = 0.007_dp+0.00005_dp*real(i,dp)
      return_value = 0.0003_dp+scale_value*(sin(0.43_dp*real(i,dp)) + &
         0.45_dp*cos(0.13_dp*real(i,dp)))
      prices(i) = prices(i-1)*exp(return_value)
   end do
   options = ufrisk_options()
   options%n_out = 10
   options%smooth = ufrisk_smooth_lpr
   options%distribution = ufrisk_distribution_student
   options%max_fit_iterations = 30
   options%smoothing_iterations = 4
   options%fractional_terms = 35
   options%log_filter_lag = 25
   options%smoothing_start = 0.18_dp

   options%model = ufrisk_model_sgarch
   short_result = varcast(prices,options)
   call assert_true(short_result%status==ufrisk_ok,'short-memory semiparametric status')
   call assert_true(allocated(short_result%short_memory_smooth%estimate),'short smoother retained')
   call assert_true(index(short_result%model_name,'Semi-')==1,'short model label')
   call assert_true(short_result%degrees_freedom>2.0_dp,'short student df')

   options%model = ufrisk_model_filgarch
   long_result = varcast(prices,options)
   call assert_true(long_result%status==ufrisk_ok,'long-memory semiparametric status')
   call assert_true(allocated(long_result%long_memory_smooth%estimate),'long smoother retained')
   call assert_true(long_result%long_memory_smooth%iterations>=1,'long smoother iterations')
   call assert_true(long_result%degrees_freedom>2.0_dp,'long student df')
   call assert_vector_finite(long_result%expected_shortfall,'long-memory ES finite')
   call finish_tests('test_semiparametric')
end program test_semiparametric
