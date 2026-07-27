! SPDX-License-Identifier: LGPL-3.0-or-later
program test_inference
   use sharper, only: dp, sr_result, sropt_result, test_result
   use sharper, only: fit_sr, fit_sropt, sr_test, sropt_test
   use sharper, only: f_ncp_krs, f_ncp_unbiased, infer_sropt, sropt_confint
   use sharper, only: power_sr_test, required_n_sr_test, power_sropt_test
   implicit none
   real(dp) :: x(80), data(80,3), ci(2), p1, p2
   integer :: i, nreq
   type(sr_result) :: z
   type(sropt_result) :: zo
   type(test_result) :: tr

   do i = 1, 80
      x(i) = 0.01_dp+0.02_dp*sin(0.37_dp*real(i,dp))+0.006_dp*cos(0.13_dp*real(i,dp))
      data(i,1) = x(i)
      data(i,2) = 0.006_dp+0.8_dp*x(i)+0.008_dp*sin(0.71_dp*real(i,dp))
      data(i,3) = 0.004_dp+0.4_dp*x(i)+0.012_dp*cos(0.43_dp*real(i,dp))
   end do
   z = fit_sr(x,ope=252.0_dp,higher_order=.true.)
   tr = sr_test(z,0.0_dp,'greater','exact')
   if (tr%p_value < 0.0_dp .or. tr%p_value > 1.0_dp) error stop 1
   if (tr%estimate <= 0.0_dp) error stop 1

   zo = fit_sropt(data,ope=252.0_dp)
   tr = sropt_test(zo,0.0_dp,'greater')
   if (tr%p_value < 0.0_dp .or. tr%p_value > 1.0_dp) error stop 1
   if (infer_sropt(zo,'krs') < 0.0_dp) error stop 1
   ci = sropt_confint(zo,0.90_dp)
   if (ci(1) > ci(2)) error stop 1

   if (f_ncp_krs(2.0_dp,4.0_dp,20.0_dp) < f_ncp_unbiased(2.0_dp,4.0_dp,20.0_dp)) error stop 1
   p1 = power_sr_test(100,0.5_dp,0.05_dp,.false.)
   p2 = power_sr_test(300,0.5_dp,0.05_dp,.false.)
   if (p2 <= p1) error stop 1
   nreq = required_n_sr_test(0.5_dp,0.8_dp,0.05_dp,.false.)
   if (power_sr_test(nreq,0.5_dp,0.05_dp,.false.) < 0.8_dp) error stop 1
   if (power_sropt_test(3,200,0.8_dp,0.05_dp) <= &
       power_sropt_test(3,60,0.8_dp,0.05_dp)) error stop 1

   print '(a)', 'test_inference: PASS'
end program test_inference
