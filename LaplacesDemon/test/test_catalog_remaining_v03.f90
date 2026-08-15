module catalog_remaining_target
use laplacesdemon, only: dp
implicit none
contains
function target(x) result(lp)
   real(dp), intent(in) :: x(:)
   real(dp) :: lp
   lp=-0.5_dp*((x(1)-1.2_dp)**2/0.7_dp+(x(2)+0.4_dp)**2/0.5_dp)
end function target
end module catalog_remaining_target

program test_catalog_remaining_v03
use laplacesdemon
use catalog_remaining_target
implicit none
integer :: fails
real(dp) :: x0(2),m(2)
real(dp), allocatable :: fw(:)
type(optim_result_t) :: o
type(mcmc_result_t) :: r

fails=0
x0=[-1.0_dp,1.0_dp]
call levenberg_marquardt_maximize(target,x0,o)
call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<1.0e-3_dp,'LM',fails)
call rprop_maximize(target,x0,o)
call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<1.0e-3_dp,'Rprop',fails)
call sgd_maximize(target,x0,o,max_iter=5000,learning_rate=0.15_dp)
call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<1.0e-3_dp,'SGD optimizer',fails)
call seed_rng(30551)
call soma_maximize(target,x0,o,max_iter=100,pop_size=20)
call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<0.08_dp,'SOMA',fails)
call seed_rng(30552)
call hit_and_run_maximize(target,x0,o,max_iter=5000)
call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<0.03_dp,'HAR optimizer',fails)

call seed_rng(30553)
call charm_sample(target,x0,[0.7_dp,0.7_dp],2400,400,2,r)
call chain_mean(r%chain,m)
call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.18_dp,'CHARM',fails)
call seed_rng(30554)
call adaptive_griddy_gibbs_sample(target,x0,[2.0_dp,2.0_dp],41,1600,300,2,r,50,fw)
call chain_mean(r%chain,m)
call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.18_dp,'AGG',fails)
call check(all(fw>0.0_dp),'AGG adaptation',fails)
call seed_rng(30555)
call sgld_sample(target,x0,0.08_dp,4000,500,5,r,0.55_dp)
call check(all(r%chain==r%chain) .and. r%acceptance_rate==1.0_dp,'SGLD',fails)
call seed_rng(30556)
call rss_sample(target,x0,[1.5_dp,1.5_dp],[-5.0_dp,-5.0_dp], &
   [5.0_dp,5.0_dp],1800,300,2,r)
call chain_mean(r%chain,m)
call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.18_dp,'RSS',fails)

if(fails==0) then
   print '(a)','test_catalog_remaining_v03: PASS'
else
   print '(a,i0)','test_catalog_remaining_v03: FAIL ',fails
   error stop 1
end if
contains
subroutine check(ok,name,nfail)
   logical,intent(in)::ok
   character(*),intent(in)::name
   integer,intent(inout)::nfail
   if(.not.ok) then
      print '(a,a)','FAIL: ',trim(name)
      nfail=nfail+1
   end if
end subroutine check
subroutine chain_mean(x,meanv)
   real(dp),intent(in)::x(:,:)
   real(dp),intent(out)::meanv(:)
   meanv=sum(x,dim=1)/real(size(x,1),dp)
end subroutine chain_mean
end program test_catalog_remaining_v03
