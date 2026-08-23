program test_parity_v03
  use mev_kinds, only: dp
  use mev_distributions, only: qgev, qgp
  use mev_bias, only: gev_bias, gev_bcor, bias_correction_result
  use mev_threshold, only: thselect_bab, bab_result
  use mev_tem, only: tem_profile_result, gpd_tem_profile, gev_tem_profile
  use mev_sampling, only: rexstudspec, rbrspec, rmev, rparp, rgparp
  use mev_mgp_likelihood, only: mgp_lik_result, gpd_to_pareto_matrix, &
    mgp_ll_log, mgp_ll_neglog, mgp_ll_br, mgp_ll_xstud, &
    mgp_cll_log, mgp_cll_neglog, mgp_cll_br, mgp_cll_xstud
  implicit none
  integer :: fails,i,n,info
  real(dp) :: b1(3),b2(3),b3(3)
  real(dp),allocatable :: x(:),xg(:),dat(:,:),tp(:,:),sp(:,:),mx(:,:),rp(:,:)
  real(dp) :: loc3(3),sc3(3),sh3(3),lam3(3),lmat(3,3),sig3(3,3),mth3(3),ll0,ar
  logical :: ok
  type(bias_correction_result) :: bc
  type(bab_result) :: br
  type(tem_profile_result) :: tr
  type(mgp_lik_result) :: ml,ml2

  fails=0
  call gev_bias([0.0_dp,1.0_dp,0.10_dp],100,b1,72)
  call gev_bias([5.0_dp,2.0_dp,0.10_dp],100,b2,72)
  call gev_bias([0.0_dp,1.0_dp,0.10_dp],200,b3,72)
  call check(all(abs(b1)<1.0_dp),'GEV Cox-Snell bias finite',fails)
  call check(abs(b2(1)-2.0_dp*b1(1))<2.0e-4_dp.and.abs(b2(2)-2.0_dp*b1(2))<2.0e-4_dp, &
    'GEV bias location-scale equivariance',fails)
  call check(abs(b2(3)-b1(3))<2.0e-4_dp,'GEV shape bias scale invariance',fails)
  call check(maxval(abs(b3-0.5_dp*b1))<2.0e-7_dp,'GEV bias inverse sample-size scaling',fails)

  n=80;allocate(xg(n),x(n))
  do i=1,n
    xg(i)=qgev((real(i,dp)-0.5_dp)/real(n,dp),loc=0.4_dp,scale=1.1_dp,shape=0.08_dp)
  end do
  call gev_bcor([0.4_dp,1.1_dp,0.08_dp],xg,bc,'subtract',nquad=48)
  call check(all(abs(bc%residual)<2.0e-3_dp),'GEV implicit bias correction solve',fails)

  do i=1,n
    x(i)=1.0_dp+qgp((real(i,dp)-0.5_dp)/real(n,dp),scale=1.0_dp,shape=0.25_dp)
  end do
  call thselect_bab(x,br,kmin=15,kmax=50,rho=-1.0_dp,test=.true.,nsim=39,level=0.90_dp)
  call check(br%convergence==0.and.allocated(br%stat),'BAB Monte Carlo test branch',fails)
  call check(size(br%stat)==br%k0_lth-1.and.all(br%ci_upper>=br%ci_lower), &
    'BAB simultaneous envelope dimensions',fails)
  call check(br%test_size>=0.0_dp.and.br%test_size<=1.0_dp,'BAB test size valid',fails)

  do i=1,n
    x(i)=qgp((real(i,dp)-0.5_dp)/real(n,dp),scale=1.2_dp,shape=0.10_dp)
  end do
  call gpd_tem_profile(x,[-0.05_dp,0.10_dp,0.25_dp],'shape',tr)
  call check(tr%convergence==0.and.all(abs(tr%r)<20.0_dp),'GPD TEM likelihood roots',fails)
  call check(all(abs(tr%rstar)<100.0_dp),'GPD modified likelihood roots finite',fails)
  call gev_tem_profile(xg,[-0.05_dp,0.08_dp,0.20_dp],'shape',tr)
  call check(tr%convergence==0.and.all(abs(tr%rstar)<100.0_dp),'GEV TEM profile correction',fails)

  allocate(dat(4,3),tp(4,3))
  tp=reshape([2.0_dp,3.0_dp,4.0_dp,5.0_dp, &
              2.5_dp,4.0_dp,3.0_dp,6.0_dp, &
              3.0_dp,2.2_dp,5.0_dp,4.5_dp],[4,3])
  dat=log(tp);loc3=0.0_dp;sc3=1.0_dp;sh3=0.0_dp;lam3=1.0_dp
  call gpd_to_pareto_matrix(dat,loc3,sc3,sh3,lam3,tp,ok)
  call check(ok.and.all(tp>1.0_dp),'multivariate GPD-to-Pareto transform',fails)
  call mgp_ll_log(dat,0.0_dp,loc3,sc3,sh3,0.5_dp,ml)
  call check(ml%convergence==0.and.ml%exponent_measure>0.0_dp,'logistic MGP likelihood',fails)
  call mgp_ll_neglog(dat,0.0_dp,loc3,sc3,sh3,1.5_dp,ml)
  call check(ml%convergence==0.and.ml%exponent_measure>0.0_dp,'negative-logistic MGP likelihood',fails)

  lmat=reshape([0.0_dp,0.4_dp,0.5_dp,0.4_dp,0.0_dp,0.35_dp,0.5_dp,0.35_dp,0.0_dp],[3,3])
  call mgp_ll_br(dat,0.0_dp,loc3,sc3,sh3,lmat,ml,nqmc=2048)
  call check(ml%convergence==0.and.ml%exponent_measure>0.0_dp,'Brown-Resnick MGP likelihood',fails)
  sig3=reshape([1.0_dp,0.3_dp,0.2_dp,0.3_dp,1.0_dp,0.25_dp,0.2_dp,0.25_dp,1.0_dp],[3,3])
  call mgp_ll_xstud(dat,0.0_dp,loc3,sc3,sh3,sig3,4.0_dp,ml,nqmc=2048)
  call check(ml%convergence==0.and.ml%exponent_measure>0.0_dp,'extremal-t MGP likelihood',fails)

  mth3=-1.0_dp
  call mgp_ll_log(dat,0.0_dp,loc3,sc3,sh3,0.5_dp,ml);ll0=ml%loglik
  call mgp_cll_log(dat,0.0_dp,mth3,loc3,sc3,sh3,0.5_dp,ml2)
  call check(ml2%convergence==0.and.abs(ml2%loglik-ll0)<1.0e-10_dp, &
    'censored logistic reduces to uncensored likelihood',fails)
  call mgp_ll_neglog(dat,0.0_dp,loc3,sc3,sh3,1.5_dp,ml);ll0=ml%loglik
  call mgp_cll_neglog(dat,0.0_dp,mth3,loc3,sc3,sh3,1.5_dp,ml2)
  call check(ml2%convergence==0.and.abs(ml2%loglik-ll0)<1.0e-9_dp, &
    'censored negative-logistic reduces to uncensored likelihood',fails)
  call mgp_ll_br(dat,0.0_dp,loc3,sc3,sh3,lmat,ml,nqmc=2048);ll0=ml%loglik
  call mgp_cll_br(dat,0.0_dp,mth3,loc3,sc3,sh3,lmat,ml2,nqmc=2048)
  call check(ml2%convergence==0.and.abs(ml2%loglik-ll0)<2.0e-8_dp, &
    'censored Brown-Resnick reduces to uncensored likelihood',fails)
  call mgp_ll_xstud(dat,0.0_dp,loc3,sc3,sh3,sig3,4.0_dp,ml,nqmc=2048);ll0=ml%loglik
  call mgp_cll_xstud(dat,0.0_dp,mth3,loc3,sc3,sh3,sig3,4.0_dp,ml2,nqmc=2048)
  call check(ml2%convergence==0.and.abs(ml2%loglik-ll0)<2.0e-8_dp, &
    'censored extremal-t reduces to uncensored likelihood',fails)

  mth3=log(3.0_dp)
  call mgp_cll_log(dat,log(3.0_dp),mth3,loc3,sc3,sh3,0.5_dp,ml)
  call check(ml%convergence==0,'censored logistic MGP likelihood',fails)
  call mgp_cll_neglog(dat,log(3.0_dp),mth3,loc3,sc3,sh3,1.5_dp,ml)
  call check(ml%convergence==0,'censored negative-logistic MGP likelihood',fails)
  call mgp_cll_br(dat,log(3.0_dp),mth3,loc3,sc3,sh3,lmat,ml,nqmc=4096)
  call check(ml%convergence==0,'censored Brown-Resnick MGP likelihood',fails)
  call mgp_cll_xstud(dat,log(3.0_dp),mth3,loc3,sc3,sh3,sig3,4.0_dp,ml,nqmc=4096)
  call check(ml%convergence==0,'censored extremal-t MGP likelihood',fails)
  call mgp_cll_log(dat,log(3.0_dp),mth3,loc3,sc3,sh3,0.5_dp,ml,likt='pois',ntot=20)
  call check(ml%convergence==0,'censored Poisson likelihood contribution',fails)


  allocate(sp(80,3),mx(120,3),rp(60,3))
  call rexstudspec(80,sig3,4.0_dp,sp,info)
  call check(info==0.and.maxval(abs(sum(sp,dim=2)-1.0_dp))<1.0e-12_dp.and.all(sp>=0.0_dp), &
    'extremal-t spectral sampler',fails)
  call rbrspec(80,sig3,sp,info)
  call check(info==0.and.maxval(abs(sum(sp,dim=2)-1.0_dp))<1.0e-12_dp.and.all(sp>0.0_dp), &
    'Brown-Resnick spectral sampler',fails)
  call rmev(120,'xstud',[4.0_dp],mx,info=info,sigma=sig3)
  call check(info==0.and.all(mx>0.0_dp),'extremal-t max-stable simulation',fails)
  call rmev(120,'br',sample=mx,info=info,sigma=sig3)
  call check(info==0.and.all(mx>0.0_dp),'Brown-Resnick max-stable simulation',fails)
  call rparp(60,1.0_dp,'sum','log',rp,par=[2.0_dp],accept_rate=ar,info=info)
  call check(info==0.and.minval(sum(rp,dim=2))>1.0_dp.and.ar>0.99_dp,'R-Pareto sum simulation',fails)
  call rgparp(60,[0.0_dp],0.0_dp,'max','log',loc3,sc3,rp,par=[2.0_dp],accept_rate=ar,info=info)
  call check(info==0.and.all(maxval(rp,dim=2)>0.0_dp).and.ar>0.0_dp,'generalized R-Pareto rejection simulation',fails)

  if(fails==0)then
    print '(a)','test_parity_v03: PASS'
  else
    print '(a,i0)','test_parity_v03: FAIL ',fails
    error stop 1
  end if
contains
  subroutine check(test,name,nfail)
    logical,intent(in)::test
    character(len=*),intent(in)::name
    integer,intent(inout)::nfail
    if(.not.test)then
      print '(a,a)','FAIL: ',trim(name);nfail=nfail+1
    end if
  end subroutine check
end program test_parity_v03
