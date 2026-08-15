module example_target
use laplacesdemon, only: dp
implicit none
contains
function log_posterior(theta) result(lp)
   real(dp),intent(in)::theta(:)
   real(dp)::lp
   lp=-0.5_dp*((theta(1)-2.0_dp)**2/1.5_dp+(theta(2)+0.5_dp)**2/0.7_dp)
end function log_posterior
end module example_target

program basic
use laplacesdemon
use example_target
implicit none
type(laplace_result_t) :: la
type(mcmc_result_t) :: mc
real(dp) :: pcov(2,2),mean(2),cov(2,2)
call seed_rng(2026)
call laplace_approximation(log_posterior,[0.0_dp,0.0_dp],la)
pcov=reshape([0.7_dp,0.0_dp,0.0_dp,0.5_dp],[2,2])
call adaptive_metropolis_sample(log_posterior,[0.0_dp,0.0_dp],pcov,7000,1000,2,mc)
call sample_covariance(mc%chain,cov,mean)
print '(a,2f10.4)', 'Laplace mode: ',la%mode
print '(a,2f10.4)', 'MCMC mean:    ',mean
print '(a,f10.4)', 'Acceptance:   ',mc%acceptance_rate
end program basic
