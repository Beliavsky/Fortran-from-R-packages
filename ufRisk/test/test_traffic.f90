! SPDX-License-Identifier: GPL-3.0-only
program test_traffic
   use ufrisk
   use test_support
   implicit none
   type(ufrisk_result) :: result
   type(ufrisk_traffic_result) :: traffic
   result%distribution = ufrisk_distribution_normal
   result%mean_return = 0.0_dp
   result%var_tail_probability = 0.05_dp
   result%es_tail_probability = 0.10_dp
   allocate(result%returns_out(10),result%sigma_forecast(10),result%var_es_level(10), &
      result%var_var_level(10),result%expected_shortfall(10))
   result%returns_out = [-0.01_dp,-0.03_dp,0.01_dp,-0.025_dp,-0.06_dp,0.0_dp, &
      -0.04_dp,-0.08_dp,-0.015_dp,-0.05_dp]
   result%sigma_forecast = 0.02_dp
   result%var_es_level = 0.03_dp
   result%var_var_level = 0.04_dp
   result%expected_shortfall = 0.05_dp
   traffic = trafftest(result)
   call assert_true(traffic%status==ufrisk_ok,'traffic status')
   call assert_true(traffic%violations_es_var==4,'ES-VaR violations')
   call assert_true(traffic%violations_var==3,'VaR violations')
   call assert_true(traffic%violations_es==2,'ES violations')
   call assert_close(traffic%breach_sum,3.696586334525814_dp,2.0e-12_dp,'breach sum')
   call assert_close(traffic%p_es_var,0.9983650626_dp,2.0e-12_dp,'traffic p ES-VaR')
   call assert_close(traffic%p_var,0.9989715020621094_dp,2.0e-12_dp,'traffic p VaR')
   call assert_close(traffic%p_es,0.9999999957121685_dp,2.0e-11_dp,'traffic p ES')
   call assert_close(traffic%weighted_absolute_deviation,14.393172669051628_dp,2.0e-12_dp,'WAD')
   call finish_tests('test_traffic')
end program test_traffic
