module core_targets
use laplacesdemon, only: dp
implicit none
contains
function quad_target(theta) result(lp)
   real(dp),intent(in)::theta(:)
   real(dp)::lp
   lp=-0.5_dp*((theta(1)-1.0_dp)**2/2.0_dp+(theta(2)+2.0_dp)**2/0.5_dp)
end function quad_target
end module core_targets

program test_core
use laplacesdemon
use core_targets
implicit none
integer :: fails,info
real(dp) :: x2(2),h(2,2),nodes(8),weights(8),bb(5,40),ll(3,4)
type(laplace_result_t) :: lap
type(waic_result_t) :: wr
type(iterative_quad_result_t) :: iq
fails=0
call seed_rng(12345)
call check(abs(normal_cdf(normal_quantile(0.975_dp,0.0_dp,1.0_dp),0.0_dp,1.0_dp)-0.975_dp)<2e-7_dp,&
   'normal cdf/quantile',fails)
call check(abs(plaplace(qlaplace(0.2_dp,1.0_dp,2.0_dp),1.0_dp,2.0_dp)-0.2_dp)<1e-12_dp,'Laplace round trip',fails)
x2=0.0_dp
call check(abs(dmvn(x2,x2,reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2]),.true.)+log(2.0_dp*pi))<1e-12_dp,&
   'MVN log density',fails)
call numerical_hessian(quad_target,[1.0_dp,-2.0_dp],h)
call check(maxval(abs(h-reshape([-0.5_dp,0.0_dp,0.0_dp,-2.0_dp],[2,2])))<2e-5_dp,'numerical Hessian',fails)
call laplace_approximation(quad_target,[0.0_dp,0.0_dp],lap,max_iter=200,tol=1e-8_dp)
call check(maxval(abs(lap%mode-[1.0_dp,-2.0_dp]))<2e-4_dp,'Laplace mode',fails)
call check(abs(lap%covariance(1,1)-2.0_dp)<2e-3_dp .and. abs(lap%covariance(2,2)-0.5_dp)<2e-3_dp,&
   'Laplace covariance',fails)
call iterative_gauss_hermite(quad_target,[0.8_dp,-1.8_dp], &
   reshape([1.5_dp,0.0_dp,0.0_dp,0.7_dp],[2,2]),8,10,1e-8_dp,iq)
call check(maxval(abs(iq%mean-[1.0_dp,-2.0_dp]))<2e-5_dp,'iterative GH mean',fails)
call check(maxval(abs(iq%covariance-reshape([2.0_dp,0.0_dp,0.0_dp,0.5_dp],[2,2])))<2e-4_dp, &
   'iterative GH covariance',fails)
call check(abs(iq%log_normalizer-log(2.0_dp*pi))<2e-4_dp,'iterative GH normalizer',fails)
call gauss_hermite_rule(8,nodes,weights)
call check(abs(sum(weights)-sqrt(pi))<1e-12_dp,'Gauss-Hermite normalization',fails)
call check(abs(sum(weights*nodes**2)-0.5_dp*sqrt(pi))<2e-11_dp,'Gauss-Hermite second moment',fails)
ll=reshape([-1.0_dp,-1.2_dp,-0.8_dp,-1.1_dp,-0.4_dp,-0.5_dp,-0.6_dp,-0.45_dp,-2.0_dp,-1.9_dp,-2.1_dp,-2.0_dp],[3,4])
call waic(ll,wr)
call check(wr%waic>0.0_dp .and. wr%p_waic>=0.0_dp,'WAIC finite',fails)
call bayesian_bootstrap_weights(5,40,bb)
call check(maxval(abs(sum(bb,dim=1)-1.0_dp))<1e-12_dp .and. minval(bb)>=0.0_dp,'Bayesian bootstrap weights',fails)
call check(abs(dpareto(2.0_dp,3.0_dp)-3.0_dp/16.0_dp)<1e-14_dp,'Pareto density',fails)
call check(abs(ddirichlet([0.2_dp,0.3_dp,0.5_dp],[1.0_dp,1.0_dp,1.0_dp])-2.0_dp)<1e-12_dp,'Dirichlet density',fails)
if(fails==0) then
   print '(a)', 'test_core: PASS'
else
   print '(a,i0)', 'test_core: FAIL count=',fails
   error stop 1
end if
contains
subroutine check(ok,name,fails)
   logical,intent(in)::ok
   character(*),intent(in)::name
   integer,intent(inout)::fails
   if(.not.ok) then; print '(a,a)', 'FAIL: ',trim(name); fails=fails+1; end if
end subroutine check
end program test_core
