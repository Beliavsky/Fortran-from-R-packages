module mcmc_targets
use laplacesdemon, only: dp
implicit none
contains
function target(theta) result(lp)
   real(dp),intent(in)::theta(:)
   real(dp)::lp
   lp=-0.5_dp*((theta(1)-1.0_dp)**2+(theta(2)+1.0_dp)**2/0.5_dp)
end function target
end module mcmc_targets

program test_mcmc
use laplacesdemon
use mcmc_targets
implicit none
integer :: fails
real(dp) :: cov(2,2),m(2),c(2,2),sir(2500,2),amcov(2,2)
type(mcmc_result_t) :: rw,am,hm,sl
fails=0
call seed_rng(777)
cov=reshape([0.7_dp,0.0_dp,0.0_dp,0.45_dp],[2,2])
call rwm_sample(target,[0.0_dp,0.0_dp],cov,12000,2000,2,rw)
call summarize(rw%chain,m,c)
call check(maxval(abs(m-[1.0_dp,-1.0_dp]))<0.10_dp,'RWM mean',fails)
call check(rw%acceptance_rate>0.15_dp .and. rw%acceptance_rate<0.75_dp,'RWM acceptance',fails)
amcov=reshape([0.2_dp,0.0_dp,0.0_dp,0.2_dp],[2,2])
call adaptive_metropolis_sample(target,[0.0_dp,0.0_dp],amcov,10000,2000,2,am)
call summarize(am%chain,m,c)
call check(maxval(abs(m-[1.0_dp,-1.0_dp]))<0.12_dp,'AM mean',fails)
call hmc_sample(target,[0.0_dp,0.0_dp],0.18_dp,8,4500,500,2,hm)
call summarize(hm%chain,m,c)
call check(maxval(abs(m-[1.0_dp,-1.0_dp]))<0.10_dp,'HMC mean',fails)
call slice_sample(target,[0.0_dp,0.0_dp],[1.0_dp,0.8_dp],4500,500,2,sl)
call summarize(sl%chain,m,c)
call check(maxval(abs(m-[1.0_dp,-1.0_dp]))<0.10_dp,'slice mean',fails)
call sir_normal(target,[0.0_dp,0.0_dp],reshape([4.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2]),8000,sir)
call summarize(sir,m,c)
call check(maxval(abs(m-[1.0_dp,-1.0_dp]))<0.14_dp,'SIR mean',fails)
call check(ess(rw%chain(:,1))>200.0_dp,'ESS positive/useful',fails)
call check(mcse_imps(rw%chain(:,1))<0.12_dp,'MCSE finite',fails)
if(fails==0) then
   print '(a)', 'test_mcmc: PASS'
else
   print '(a,i0)', 'test_mcmc: FAIL count=',fails
   error stop 1
end if
contains
subroutine summarize(x,m,c)
   real(dp),intent(in)::x(:,:)
   real(dp),intent(out)::m(:),c(:,:)
   call sample_covariance(x,c,m)
end subroutine summarize
subroutine check(ok,name,fails)
   logical,intent(in)::ok
   character(*),intent(in)::name
   integer,intent(inout)::fails
   if(.not.ok) then; print '(a,a)', 'FAIL: ',trim(name); fails=fails+1; end if
end subroutine check
end program test_mcmc
