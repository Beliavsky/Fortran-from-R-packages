module advanced_example_target
use laplacesdemon, only: dp
implicit none
contains
function log_target(theta) result(lp)
   real(dp), intent(in) :: theta(:)
   real(dp) :: lp
   lp=-0.5_dp*((theta(1)-1.5_dp)**2/0.8_dp+(theta(2)+0.4_dp)**2/0.5_dp)
end function log_target
end module advanced_example_target

program advanced_v02
use laplacesdemon
use advanced_example_target
implicit none
type(mcmc_result_t) :: nuts
type(vb_result_t) :: vb
real(dp) :: mean_chain(2),step
call seed_rng(24680)
call nuts_sample(log_target,[0.0_dp,0.0_dp],1800,600,2,nuts,adapt_steps=600,final_step=step)
mean_chain=sum(nuts%chain,dim=1)/real(size(nuts%chain,1),dp)
call seed_rng(24681)
call variational_bayes_salimans2(log_target,[0.0_dp,0.0_dp],vb,iterations=700)
print '(a,2f10.4)', 'NUTS posterior mean: ',mean_chain
print '(a,f10.4)', 'NUTS final step size: ',step
print '(a,2f10.4)', 'VB mean: ',vb%mean
end program advanced_v02
