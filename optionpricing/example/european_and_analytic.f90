! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program european_and_analytic
   use optionpricing, only : dp, european_result, moments_result, bs_ec, bs_ep, &
      eval_ecv, eval_lb, eval_eqcv, asian_call_app_lord
   implicit none
   type(european_result) :: call,put
   type(moments_result) :: ecv,lb,qcv

   call=bs_ec(0.25_dp,100.0_dp,0.05_dp,0.2_dp,100.0_dp)
   put=bs_ep(0.25_dp,100.0_dp,0.05_dp,0.2_dp,100.0_dp)
   ecv=eval_ecv(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp)
   lb=eval_lb(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp,.false.)
   qcv=eval_eqcv(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp)

   print '(a,3(f12.6,1x))','European call: ',call%price,call%delta,call%gamma
   print '(a,3(f12.6,1x))','European put:  ',put%price,put%delta,put%gamma
   print '(a,3(f12.6,1x))','Curran ECV:    ',ecv%price,ecv%delta,ecv%gamma
   print '(a,3(es12.4,1x))','Lower-bound CV:',lb%price,lb%delta,lb%gamma
   print '(a,3(es12.4,1x))','Quadratic CV:  ',qcv%price,qcv%delta,qcv%gamma
   print '(a,f12.6)','Lord approximation: ', &
      asian_call_app_lord(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp,.true.)
end program european_and_analytic
