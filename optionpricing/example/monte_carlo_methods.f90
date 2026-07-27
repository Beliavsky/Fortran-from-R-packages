! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program monte_carlo_methods
   use optionpricing, only : dp, greeks_result, asian_call_naive_mc, &
      asian_call_best_mc, asian_call_naive_qmc, asian_call_best_qmc
   implicit none
   type(greeks_result) :: naive_mc,best_mc,naive_qmc,best_qmc

   naive_mc=asian_call_naive_mc(5000,1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp,123)
   best_mc=asian_call_best_mc(5000,1.0_dp,12,100.0_dp,0.05_dp,0.2_dp, &
      100.0_dp,500,80,1.0e-13_dp,123)
   naive_qmc=asian_call_naive_qmc(16,257,76,1.0_dp,12,100.0_dp,0.05_dp, &
      0.2_dp,100.0_dp,'pca',.true.,123)
   best_qmc=asian_call_best_qmc(16,257,76,1.0_dp,12,100.0_dp,0.05_dp, &
      0.2_dp,100.0_dp,'pca',1,.true.,'splitting',123,maxiter=80,tol=1.0e-13_dp)

   call show('Naive MC ',naive_mc)
   call show('Best MC  ',best_mc)
   call show('Naive QMC',naive_qmc)
   call show('Best QMC ',best_qmc)
contains
   subroutine show(label,res)
      character(len=*), intent(in) :: label
      type(greeks_result), intent(in) :: res
      print '(a,1x,3(f12.6,1x))',label,res%estimate
      print '(a,1x,3(es12.4,1x))','95% errors',res%error95
   end subroutine show
end program monte_carlo_methods
