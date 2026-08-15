module completion_target
use laplacesdemon, only: dp
implicit none
contains
function log_target(x) result(lp)
   real(dp), intent(in) :: x(:)
   real(dp) :: lp
   lp=-0.5_dp*((x(1)-1.2_dp)**2/0.7_dp+(x(2)+0.4_dp)**2/0.5_dp)
end function log_target
end module completion_target

program completion_v03
use laplacesdemon
use completion_target
implicit none
type(optim_result_t) :: opt
type(iterative_quad_result_t) :: quad
type(mcmc_result_t) :: chain
real(dp) :: x0(2),cov(2,2),post_mean(2)

x0=[-1.0_dp,1.0_dp]
cov=0.0_dp
cov(1,1)=1.0_dp
cov(2,2)=1.0_dp

call trust_region_maximize(log_target,x0,opt)
call adaptive_sparse_grid_quadrature(log_target,[0.0_dp,0.0_dp],cov,4,8,1.0e-6_dp,quad)
call seed_rng(30399)
call admg_sample(log_target,x0,cov,1600,300,2,chain,periodicity=25)
post_mean=sum(chain%chain,dim=1)/real(size(chain%chain,1),dp)

print '(a,2f10.4)', 'Trust-region mode:      ',opt%par
print '(a,2f10.4)', 'Sparse-grid mean:       ',quad%mean
print '(a,2f10.4)', 'ADMG posterior mean:    ',post_mean
print '(a,f10.4)',  'ADMG acceptance rate:  ',real(chain%accepted,dp)/real(chain%proposed,dp)
end program completion_v03
