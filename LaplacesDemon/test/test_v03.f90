module v03_targets
use laplacesdemon, only: dp
implicit none
contains
function target(x) result(lp)
   real(dp), intent(in) :: x(:)
   real(dp) :: lp
   lp=-0.5_dp*((x(1)-1.2_dp)**2/0.7_dp+(x(2)+0.4_dp)**2/0.5_dp)
end function target
function target_positive(x) result(lp)
   real(dp), intent(in) :: x(:)
   real(dp) :: lp
   if(any(x<=0.0_dp)) then; lp=-huge(1.0_dp); else; lp=-0.5_dp*sum((x-[1.2_dp,0.8_dp])**2/0.2_dp); end if
end function target_positive
subroutine components(x,y)
   real(dp), intent(in) :: x(:)
   real(dp), intent(out) :: y(:)
   integer :: i
   do i=1,size(y); y(i)=-0.5_dp*(x(1)-(0.5_dp+0.02_dp*real(i,dp)))**2; end do
end subroutine components
end module v03_targets

program test_v03
use laplacesdemon
use v03_targets
implicit none
integer :: fails,i,j
real(dp) :: x0(2),cov(2,2),m(2),v(2),temps(3),swap
real(dp) :: walkers(3,2),stats(2),pvals(2),wstick(5),prec(2,2)
real(dp), allocatable :: coupled(:,:,:),hd(:,:),fw(:)
type(optim_result_t) :: o
type(mcmc_result_t) :: r
type(iterative_quad_result_t) :: iq
type(heidelberger_result_t) :: hw
real(dp), allocatable :: diagchain(:,:)
fails=0; x0=[-1.0_dp,1.0_dp]; cov=0.0_dp; cov(1,1)=1.0_dp; cov(2,2)=1.0_dp

call newton_maximize(target,x0,o); call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<2.0e-4_dp,'Newton',fails)
call nelder_mead_maximize(target,x0,o,max_iter=500); call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<2.0e-3_dp,'NM',fails)
call conjugate_gradient_maximize(target,x0,o); call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<2.0e-3_dp,'CG',fails)
call dfp_maximize(target,x0,o); call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<2.0e-3_dp,'DFP',fails)
call lbfgs_maximize(target,x0,o); call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<2.0e-3_dp,'LBFGS',fails)
call hooke_jeeves_maximize(target,x0,o); call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<5.0e-3_dp,'HJ',fails)
call trust_region_maximize(target,x0,o); call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<5.0e-3_dp,'TR',fails)
call spg_maximize(target,x0,o); call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<5.0e-3_dp,'SPG',fails)
call sr1_maximize(target,x0,o); call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<5.0e-3_dp,'SR1',fails)
call bhhh_maximize(components,[0.0_dp],50,o); call check(abs(o%par(1)-1.01_dp)<0.03_dp,'BHHH',fails)
call seed_rng(30301); call pso_maximize(target,x0,o,max_iter=250,n_particles=20)
call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<0.05_dp,'PSO',fails)
call seed_rng(30302); call genetic_maximize(target,x0,o,max_iter=250,pop_size=30)
call check(maxval(abs(o%par-[1.2_dp,-0.4_dp]))<0.08_dp,'AGA',fails)

call seed_rng(30310); call independence_metropolis_sample(target,x0,[1.0_dp,-0.2_dp],cov,2500,500,2,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.15_dp,'IM',fails)
call seed_rng(30311); call multiple_try_metropolis_sample(target,x0,0.6_dp*cov,4,2200,400,2,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.18_dp,'MTM',fails)
call seed_rng(30312); call hit_and_run_sample(target,x0,0.9_dp,2500,500,2,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.18_dp,'HARM',fails)
call seed_rng(30313); call adaptive_mwg_sample(target,x0,[0.5_dp,0.5_dp],2500,500,2,r,final_sd=fw)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.18_dp,'AMWG',fails)
call seed_rng(30314); call adaptive_mixture_metropolis_sample(target,x0,cov,2500,500,2,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.18_dp,'AMM',fails)
call seed_rng(30315); call griddy_gibbs_sample(target,x0,[2.0_dp,2.0_dp],41,1400,300,2,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.15_dp,'GG',fails)
call seed_rng(30316); call ohss_sample(target,x0,1.5_dp,1800,300,2,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.18_dp,'OHSS',fails)
call seed_rng(30317); call uess_sample(target,x0,cov,1.5_dp,1800,300,2,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.18_dp,'UESS',fails)
call seed_rng(30318); call afss_sample(target,x0,cov,1.5_dp,1800,300,2,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.18_dp,'AFSS',fails)
call seed_rng(30319); call tempered_hmc_sample(target,x0,0.15_dp,8,1.5_dp,1800,300,2,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,-0.4_dp]))<0.20_dp,'THMC',fails)
call seed_rng(30320); call random_dive_sample(target_positive,[1.0_dp,1.0_dp],0.25_dp,2500,500,2,r)
call moments(r%chain,m,v); call check(maxval(abs(m-[1.2_dp,0.8_dp]))<0.18_dp,'RDMH',fails)
walkers=reshape([0.0_dp,1.0_dp,2.0_dp,0.0_dp,-0.5_dp,-1.0_dp],[3,2]); temps=[1.0_dp,2.0_dp,4.0_dp]
call seed_rng(30321); call mcmcmc_sample(target,walkers,temps,1200,200,2,coupled,swap,0.5_dp)
call check(swap>=0.0_dp .and. swap<=1.0_dp,'MCMCMC swaps',fails)

call componentwise_iterative_quadrature(target,[0.0_dp,0.0_dp],[1.0_dp,1.0_dp],9,20,1.0e-8_dp,iq)
call check(maxval(abs(iq%mean-[1.2_dp,-0.4_dp]))<2.0e-3_dp,'CAGH mean',fails)
call adaptive_sparse_grid_quadrature(target,[0.0_dp,0.0_dp],cov,4,8,1.0e-6_dp,iq)
call check(maxval(abs(iq%mean-[1.2_dp,-0.4_dp]))<0.03_dp,'AGHSG mean',fails)

allocate(diagchain(1000,2)); call seed_rng(30330)
do i=1,1000; diagchain(i,1)=rand_normal(); diagchain(i,2)=0.5_dp*rand_normal(); end do
call bmk_diagnostic(diagchain,10,hd); call check(maxval(hd)<0.65_dp,'BMK',fails)
call ks_diagnostic(diagchain,stats,pvals); call check(minval(pvals)>1.0e-4_dp,'KS',fails)
call heidelberger_diagnostic(diagchain,hw); call check(any(hw%stationarity_pass),'Heidelberger',fails)

prec=0.0_dp; prec(1,1)=2.0_dp; prec(2,2)=3.0_dp
call check(abs(dmvnormal_precision([0.2_dp,-0.1_dp],[0.0_dp,0.0_dp],prec,.true.)-&
     dmvn([0.2_dp,-0.1_dp],[0.0_dp,0.0_dp],reshape([0.5_dp,0.0_dp,0.0_dp,1.0_dp/3.0_dp],[2,2]),.true.))<1.0e-10_dp,&
     'MVN precision',fails)
call check(abs(dmvc([0.2_dp,-0.1_dp],[0.0_dp,0.0_dp],cov,.true.) &
     -dmvt([0.2_dp,-0.1_dp],[0.0_dp,0.0_dp],cov,1.0_dp,.true.))<1.0e-12_dp,'MVC',fails)
call check(dmatrixgamma(reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2]),3.0_dp,1.0_dp,cov)>0.0_dp,'matrix gamma',fails)
call check(dnormlaplace(0.2_dp,0.0_dp,1.0_dp,1.2_dp,0.8_dp)>0.0_dp,'normal Laplace',fails)
call rstick(4,2.0_dp,wstick); call check(abs(sum(wstick)-1.0_dp)<1.0e-12_dp,'stick RNG',fails)

if(fails==0) then
   print '(a)', 'test_v03: PASS'
else
   print '(a,i0)', 'test_v03: FAIL ',fails
   error stop 1
end if
contains
subroutine check(ok,name,fails)
   logical,intent(in)::ok
   character(*),intent(in)::name
   integer,intent(inout)::fails
   if(.not.ok) then; print '(a,a)', 'FAIL: ',trim(name); fails=fails+1; end if
end subroutine check
subroutine moments(x,m,v)
   real(dp), intent(in) :: x(:,:)
   real(dp), intent(out) :: m(:),v(:)
   integer :: ii
   m=sum(x,dim=1)/real(size(x,1),dp); v=0.0_dp
   do ii=1,size(x,1); v=v+(x(ii,:)-m)**2; end do
   v=v/real(max(1,size(x,1)-1),dp)
end subroutine moments
end program test_v03
