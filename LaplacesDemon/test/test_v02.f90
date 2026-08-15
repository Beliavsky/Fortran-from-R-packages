module v02_targets
use laplacesdemon, only: dp
implicit none
contains
function target(x) result(lp)
   real(dp), intent(in) :: x(:)
   real(dp) :: lp,d1,d2
   d1=x(1)-1.0_dp; d2=x(2)+0.5_dp
   lp=-0.5_dp*(d1*d1/0.5_dp+d2*d2/0.4_dp)
end function target
function flat_target(x) result(lp)
   real(dp), intent(in) :: x(:)
   real(dp) :: lp
   lp=0.0_dp*sum(x)
end function flat_target
end module v02_targets

program test_v02
use laplacesdemon
use v02_targets
implicit none
real(dp) :: x0(2),cov(2,2),m(2),v(2),step
real(dp) :: xp0(2),acc
real(dp), allocatable :: fcov(:,:)
type(mcmc_result_t) :: r
type(vb_result_t) :: vb
type(raftery_result_t) :: raf
integer :: fails,disc(200),df,i
real(dp) :: hstat
fails=0
x0=[-1.0_dp,1.5_dp]
cov=0.0_dp; cov(1,1)=0.5_dp; cov(2,2)=0.5_dp

call seed_rng(20260814)
call nuts_sample(target,x0,1800,600,2,r,adapt_steps=600,target_accept=0.70_dp,max_depth=7,final_step=step)
call moments(r%chain,m,v)
call check(maxval(abs(m-[1.0_dp,-0.5_dp]))<=0.25_dp,'NUTS posterior mean',fails)
call check(maxval(abs(v-[0.5_dp,0.4_dp]))<=0.22_dp,'NUTS posterior variance',fails)
call check(step>0.0_dp .and. step<5.0_dp,'NUTS final step size',fails)

call seed_rng(20260815)
call hmcda_sample(target,x0,1600,500,2,1.2_dp,r,adapt_steps=500,target_accept=0.65_dp,final_step=step)
call moments(r%chain,m,v)
call check(maxval(abs(m-[1.0_dp,-0.5_dp]))<=0.25_dp,'HMCDA posterior mean',fails)
call check(r%acceptance_rate>=0.25_dp,'HMCDA acceptance rate',fails)

call seed_rng(20260816)
call dram_sample(target,x0,cov,2200,600,2,r,adapt_start=100,periodicity=25,final_cov=fcov)
call moments(r%chain,m,v)
call check(maxval(abs(m-[1.0_dp,-0.5_dp]))<=0.25_dp,'DRAM posterior mean',fails)
call check(fcov(1,1)>0.0_dp .and. fcov(2,2)>0.0_dp,'DRAM covariance',fails)

call seed_rng(20260817)
call ram_sample(target,x0,cov,2200,600,2,r,target_accept=0.234_dp,final_cov=fcov)
call moments(r%chain,m,v)
call check(maxval(abs(m-[1.0_dp,-0.5_dp]))<=0.28_dp,'RAM posterior mean',fails)
call check(r%acceptance_rate>=0.05_dp .and. r%acceptance_rate<=0.8_dp,'RAM acceptance rate',fails)

call seed_rng(20260818)
xp0=[2.0_dp,-1.5_dp]
call twalk_sample(target,x0,xp0,3000,800,2,r)
call moments(r%chain,m,v)
call check(maxval(abs(m-[1.0_dp,-0.5_dp]))<=0.35_dp,'t-walk posterior mean',fails)
call check(r%acceptance_rate>0.01_dp,'t-walk acceptance rate',fails)

call seed_rng(20260819)
call pcn_sample(flat_target,[0.0_dp,0.0_dp],identity(2),0.35_dp,4000,1000,3,r)
call moments(r%chain,m,v)
call check(maxval(abs(m))<=0.35_dp,'pCN prior mean',fails)
call check(all(v>=0.55_dp) .and. all(v<=1.45_dp),'pCN prior variance',fails)

call seed_rng(20260820)
call variational_bayes_salimans2(target,[0.0_dp,0.0_dp],vb,iterations=700,stop_tolerance=2.0e-3_dp)
call check(maxval(abs(vb%mean-[1.0_dp,-0.5_dp]))<=0.25_dp,'Salimans2 VB mean',fails)
call check(all(diagonal(vb%covariance)>=0.25_dp) .and. &
     all(diagonal(vb%covariance)<=0.85_dp),'Salimans2 VB covariance',fails)

call raftery_diagnostic(r%chain,raf,q=0.25_dp,r=0.06_dp,s=0.90_dp)
call check(raf%enough_samples,'Raftery sufficient sample flag',fails)
call check(all(raf%total>0),'Raftery total draws',fails)
do i=1,200
   disc(i)=mod(i,3)
end do
call hangartner_chisq(disc,4,hstat,df)
call check(df>0 .and. hstat>=0.0_dp,'Hangartner statistic',fails)

if(fails==0) then
   print '(a)', 'test_v02: PASS'
else
   print '(a,i0)', 'test_v02: FAIL ',fails
   error stop 1
end if
contains
subroutine check(ok,label,fails)
   logical, intent(in) :: ok
   character(len=*), intent(in) :: label
   integer, intent(inout) :: fails
   if(.not.ok) then
      fails=fails+1
      print '(a,a)', 'test_v02 failure: ',trim(label)
   end if
end subroutine check
subroutine moments(x,m,v)
   real(dp), intent(in) :: x(:,:)
   real(dp), intent(out) :: m(:),v(:)
   integer :: i
   m=sum(x,dim=1)/real(size(x,1),dp); v=0.0_dp
   do i=1,size(x,1); v=v+(x(i,:)-m)**2; end do
   v=v/real(max(1,size(x,1)-1),dp)
end subroutine moments
pure function identity(n) result(a)
   integer, intent(in) :: n
   real(dp) :: a(n,n)
   integer :: i
   a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
end function identity
pure function diagonal(a) result(d)
   real(dp), intent(in) :: a(:,:)
   real(dp) :: d(min(size(a,1),size(a,2)))
   integer :: i
   do i=1,size(d); d(i)=a(i,i); end do
end function diagonal
end program test_v02
