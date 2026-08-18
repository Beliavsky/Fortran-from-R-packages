program test_residual_break
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : set_seed
   use mcmcpack_changepoint, only : change_result, mcmc_residual_break_analysis
   implicit none
   real(dp) :: y(12),a0(2,2)
   type(change_result) :: r
   integer :: i
   y=[(-1.0_dp+0.05_dp*real(i,dp),i=1,6),(1.0_dp+0.05_dp*real(i-6,dp),i=7,12)]
   a0=0.0_dp;a0(1,1)=5.0_dp;a0(1,2)=1.0_dp;a0(2,2)=1.0_dp
   call set_seed(1717)
   r=mcmc_residual_break_analysis(y,1,0.0_dp,0.001_dp,0.1_dp,0.1_dp,a0,10,20,2)
   if(r%status/=0) error stop 'residual break status'
   if(any(shape(r%draws)/=[10,8])) error stop 'residual break shape'
   if(any(shape(r%states)/=[10,12])) error stop 'residual break states'
   print *,'test_residual_break: PASS'
end program
