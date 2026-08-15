module completion_targets
use laplacesdemon, only: dp, rand_normal
implicit none
contains
function normal_target(x) result(lp)
   real(dp), intent(in) :: x(:)
   real(dp) :: lp
   lp=-0.5_dp*((x(1)-1.2_dp)**2/0.7_dp+(x(2)+0.4_dp)**2/0.5_dp)
end function normal_target
subroutine normal_conditional(j,x)
   integer, intent(in) :: j
   real(dp), intent(inout) :: x(:)
   if(j==1) then
      x(1)=1.2_dp+sqrt(0.7_dp)*rand_normal()
   else
      x(2)=-0.4_dp+sqrt(0.5_dp)*rand_normal()
   end if
end subroutine normal_conditional
function sparse_target(x) result(lp)
   real(dp), intent(in) :: x(:)
   real(dp) :: lp
   lp=-0.5_dp*((x(1)-1.5_dp)**2/0.15_dp+x(2)**2/0.5_dp)
end function sparse_target
end module completion_targets

program test_v03_completion
use laplacesdemon
use completion_targets
implicit none
integer :: fails,i
real(dp) :: x0(2),cov(2,2),m(2),v(2),chains(3,2),acc
real(dp), allocatable :: out(:,:,:),fc(:,:),feps(:),fsd(:)
real(dp) :: baseline(100,2)
integer :: dyn(1)
logical :: selectable(2),initial_active(2)
real(dp) :: pprob(2)
type(mcmc_result_t) :: r
type(rj_selection_result_t) :: rr
fails=0; x0=[-0.5_dp,0.8_dp]; cov=0.0_dp; cov(1,1)=0.5_dp; cov(2,2)=0.5_dp

call seed_rng(30401); call admg_sample(normal_target,x0,cov,1800,300,2,r,periodicity=50,final_cov=fc)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.25_dp,'ADMG',fails)
call seed_rng(30402); call adaptive_hmc_sample(normal_target,x0,[0.08_dp,0.08_dp],3,3000,500,2,r,50,feps)
call moments(r%chain,m,v)
call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.20_dp,'AHMC',fails)
call seed_rng(30403); call drm_sample(normal_target,x0,cov,1800,300,2,r,0.5_dp)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.25_dp,'DRM',fails)
call seed_rng(30404); call gibbs_sample(normal_target,normal_conditional,x0,1000,100,1,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.10_dp,'Gibbs',fails)
chains=reshape([-0.5_dp,0.0_dp,0.7_dp,0.8_dp,-0.9_dp,-0.1_dp],[3,2])
call seed_rng(30405); call inca_sample(normal_target,chains,cov,1000,200,2,out,acc,50,100)
call check(acc>0.01_dp .and. acc<0.95_dp,'INCA acceptance',fails)
call check(all(out==out),'INCA finite',fails)
call seed_rng(30406); call refractive_sample(normal_target,x0,0.02_dp,1.2_dp,2,1200,200,2,r,adapt=.true.)
call moments(r%chain,m,v); call check(all(m==m) .and. r%acceptance_rate>0.0_dp,'Refractive',fails)
call seed_rng(30407); call smwg_sample(normal_target,x0,[0.5_dp,0.5_dp],1400,200,2,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.25_dp,'SMWG',fails)
call seed_rng(30408); call samwg_sample(normal_target,x0,[0.5_dp,0.5_dp],1400,300,2,r,50,fsd)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.25_dp,'SAMWG',fails)

do i=1,100
   baseline(i,1)=1.2_dp+sqrt(0.7_dp)*rand_normal()
   baseline(i,2)=-0.4_dp+sqrt(0.5_dp)*rand_normal()
end do
dyn=[2]
call seed_rng(30409); call usmwg_sample(normal_target,baseline,dyn,x0,[0.3_dp,0.3_dp],500,100,1,r)
call check(all(r%chain==r%chain),'USMWG',fails)
call seed_rng(30410); call usamwg_sample(normal_target,baseline,dyn,x0,[0.3_dp,0.3_dp],500,100,1,r)
call check(all(r%chain==r%chain),'USAMWG',fails)

selectable=.true.; initial_active=[.true.,.false.]; pprob=[0.6_dp,0.4_dp]
call seed_rng(30411)
call reversible_jump_selection_sample(sparse_target,[1.0_dp,0.0_dp],selectable,initial_active,2,pprob,0.5_dp,1200,200,2,rr)
call check(rr%acceptance_rate>0.0_dp .and. rr%acceptance_rate<1.0_dp,'RJ acceptance',fails)
call check(count(rr%active(:,1))>size(rr%active,1)/2,'RJ signal inclusion',fails)

if(fails==0) then
   print '(a)', 'test_v03_completion: PASS'
else
   print '(a,i0)', 'test_v03_completion: FAIL ',fails
   error stop 1
end if
contains
subroutine check(ok,name,nfail)
   logical,intent(in)::ok
   character(*),intent(in)::name
   integer,intent(inout)::nfail
   if(.not.ok) then; print '(a,a)', 'FAIL: ',trim(name); nfail=nfail+1; end if
end subroutine check
subroutine moments(x,meanv,varv)
   real(dp),intent(in)::x(:,:)
   real(dp),intent(out)::meanv(:),varv(:)
   integer::j
   meanv=sum(x,dim=1)/real(size(x,1),dp); varv=0.0_dp
   do j=1,size(x,1); varv=varv+(x(j,:)-meanv)**2; end do
   varv=varv/real(max(1,size(x,1)-1),dp)
end subroutine moments
end program test_v03_completion
