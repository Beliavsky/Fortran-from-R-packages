module algorithm_targets
use laplacesdemon, only: dp
implicit none
contains
function target2(theta) result(lp)
   real(dp),intent(in)::theta(:)
   real(dp)::lp
   lp=-0.5_dp*((theta(1)-1.0_dp)**2+(theta(2)+1.0_dp)**2/0.5_dp)
end function target2
function ell_lik(theta) result(lp)
   real(dp),intent(in)::theta(:)
   real(dp)::lp
   lp=-0.5_dp*(theta(1)-1.0_dp)**2
end function ell_lik
end module algorithm_targets

program test_algorithms
use laplacesdemon
use algorithm_targets
implicit none
integer :: fails,i,j,k
real(dp) :: m(2),c(2,2),mu(2),cov(2,2),essh(4)
real(dp) :: walkers(10,2),chains(8,2),acc
real(dp),allocatable :: ens(:,:,:),de(:,:,:),flat(:,:)
type(mcmc_result_t) :: mw,ma,es
fails=0
call seed_rng(9981)
call mwg_sample(target2,[0.0_dp,0.0_dp],[0.8_dp,0.6_dp],7000,1000,2,mw)
call sample_covariance(mw%chain,c,m)
call check(maxval(abs(m-[1.0_dp,-1.0_dp]))<0.12_dp,'MWG mean',fails)
call mala_sample(target2,[0.0_dp,0.0_dp],0.45_dp,6000,1000,2,ma)
call sample_covariance(ma%chain,c,m)
call check(maxval(abs(m-[1.0_dp,-1.0_dp]))<0.12_dp,'MALA mean',fails)
call elliptical_slice_sample(ell_lik,[0.0_dp],reshape([1.0_dp],[1,1]),5000,500,2,es)
call check(abs(sum(es%chain(:,1))/real(size(es%chain,1),dp)-0.5_dp)<0.08_dp,'elliptical slice mean',fails)
do i=1,10
   walkers(i,1)=1.0_dp+1.5_dp*rand_normal()
   walkers(i,2)=-1.0_dp+1.0_dp*rand_normal()
end do
call aies_sample(target2,walkers,3500,500,2,ens,acc)
allocate(flat(size(ens,1)*size(ens,2),2)); k=0
do i=1,size(ens,1); do j=1,size(ens,2); k=k+1; flat(k,:)=ens(i,j,:); end do; end do
call sample_covariance(flat,c,m)
call check(maxval(abs(m-[1.0_dp,-1.0_dp]))<0.14_dp,'AIES mean',fails)
deallocate(flat)
do i=1,8
   chains(i,1)=1.0_dp+1.5_dp*rand_normal()
   chains(i,2)=-1.0_dp+1.0_dp*rand_normal()
end do
call demc_sample(target2,chains,3500,500,2,de,acc)
allocate(flat(size(de,1)*size(de,2),2)); k=0
do i=1,size(de,1); do j=1,size(de,2); k=k+1; flat(k,:)=de(i,j,:); end do; end do
call sample_covariance(flat,c,m)
call check(maxval(abs(m-[1.0_dp,-1.0_dp]))<0.15_dp,'DEMC mean',fails)
mu=[0.0_dp,0.0_dp]; cov=reshape([4.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2])
call pmc_normal(target2,mu,cov,5000,4,mu,cov,essh)
call check(maxval(abs(mu-[1.0_dp,-1.0_dp]))<0.10_dp,'PMC mean',fails)
call check(minval(essh)>100.0_dp,'PMC ESS',fails)
if(fails==0) then
   print '(a)', 'test_algorithms: PASS'
else
   print '(a,i0)', 'test_algorithms: FAIL count=',fails
   error stop 1
end if
contains
subroutine check(ok,name,fails)
   logical,intent(in)::ok
   character(*),intent(in)::name
   integer,intent(inout)::fails
   if(.not.ok) then; print '(a,a)', 'FAIL: ',trim(name); fails=fails+1; end if
end subroutine check
end program test_algorithms
