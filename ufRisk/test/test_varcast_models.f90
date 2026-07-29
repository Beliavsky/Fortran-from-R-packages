! SPDX-License-Identifier: GPL-3.0-only
program test_varcast_models
   use ufrisk
   use test_support
   implicit none
   integer, parameter :: n = 110
   real(dp) :: prices(n),return_value
   integer :: i,model
   type(ufrisk_options) :: options
   type(ufrisk_result) :: result

   prices(1) = 100.0_dp
   do i = 2,n
      return_value = 0.0004_dp+0.009_dp*sin(0.29_dp*real(i,dp)) + &
         0.004_dp*cos(0.11_dp*real(i,dp))+0.002_dp*sin(0.053_dp*real(i*i,dp))
      prices(i) = prices(i-1)*exp(return_value)
   end do
   do model = ufrisk_model_sgarch,ufrisk_model_filgarch
      options = ufrisk_options()
      options%model = model
      options%n_out = 10
      options%distribution = ufrisk_distribution_normal
      options%max_fit_iterations = 35
      options%fractional_terms = 40
      options%log_filter_lag = 30
      result = varcast(prices,options)
      call assert_true(result%status==ufrisk_ok,'parametric model status')
      call assert_true(size(result%sigma_forecast)==10,'forecast size')
      call assert_vector_finite(result%sigma_forecast,'finite sigma forecast')
      call assert_vector_finite(result%var_var_level,'finite VaR')
      call assert_vector_finite(result%expected_shortfall,'finite ES')
      call assert_true(all(result%sigma_forecast>0.0_dp),'positive sigma forecast')
      call assert_true(all(result%expected_shortfall>=result%var_es_level),'ES exceeds matching VaR')
   end do
   call finish_tests('test_varcast_models')
end program test_varcast_models
