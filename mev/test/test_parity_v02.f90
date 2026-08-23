program test_parity_v02
  use mev_kinds, only: dp
  use mev_distributions, only: qgp, qgev
  use mev_tailindex, only: shape_hill
  use mev_threshold
  use mev_bias
  use mev_stein
  use mev_erm
  use mev_mgp
  use mev_taildep_kj
  use mev_profile
  use mev_math, only: normal_cdf
  implicit none
  integer :: fails,i,n
  real(dp), allocatable :: x(:),path(:),xgev(:),umat(:,:)
  real(dp) :: h,lt,tr,b(2),sig2(2,2),up2(2),z2(2),lam2(2,2),v,expect
  real(dp) :: sig3(3,3),z3(3),qmat(3,3),lvec(3)
  type(bab_result) :: br
  type(bias_correction_result) :: bc
  type(wgpd_fit_result) :: wf
  type(erm_result) :: er
  type(kjtail_result) :: kj
  type(profile_result) :: pr
  integer :: kvals(2)

  fails=0
  n=240
  allocate(x(n))
  do i=1,n
    x(i)=qgp((real(i,dp)-0.5_dp)/real(n,dp),scale=1.0_dp,shape=0.35_dp)+1.0_dp
  end do

  h=shape_hill(x,80)
  lt=shape_lthill(x,80,80)
  tr=shape_trimhill(x,80,0)
  call check(abs(lt-h)<1.0e-12_dp,'lower-trimmed Hill reduces to Hill',fails)
  call check(abs(tr-h)<1.0e-12_dp,'trimmed Hill k0=0 reduces to Hill',fails)
  call shape_lthill_path(x,80,path,1)
  call check(size(path)==80.and.abs(path(80)-h)<1.0e-12_dp,'lower-trimmed Hill path',fails)
  call check(abs(bab_fcst(-1.0_dp)-0.006954665055_dp)<2.0e-10_dp,'BAB expint constant',fails)
  call thselect_bab(x,br,20,120,-1.0_dp)
  call check(br%convergence==0.and.br%k0>0.and.br%k0_lth>=20,'BAB threshold selection',fails)

  call gpd_bias([1.0_dp,0.0_dp],100,b)
  call check(maxval(abs(b-[0.03_dp,-0.03_dp]))<1.0e-14_dp,'GPD Cox-Snell bias formula',fails)
  do i=1,n
    x(i)=qgp((real(i,dp)-0.5_dp)/real(n,dp),scale=1.2_dp,shape=0.10_dp)
  end do
  call gpd_bcor([1.2_dp,0.10_dp],x,bc,'subtract')
  call check(bc%convergence==0.and.maxval(abs(bc%residual))<2.0e-5_dp,'GPD implicit bias correction',fails)
  call gpd_bcor([1.2_dp,0.10_dp],x,bc,'firth',.true.)
  call check(maxval(abs(bc%residual))<5.0e-3_dp,'GPD Firth score solve',fails)

  call fit_wgpd(x,wf,0.0_dp,1.0_dp)
  call check(wf%convergence==0.and.wf%estimate(1)>0.0_dp,'Stein weighted GPD fit',fails)
  call check(wf%estimate(2)>-1.0_dp.and.wf%estimate(2)<2.0_dp,'Stein shape bounds',fails)

  kvals=[50,100]
  x=x+1.0_dp
  call shape_erm(x,kvals,er,'bdgm')
  call check(all(er%shape>0.0_dp).and.all(er%rho<0.0_dp),'BDGM exponential regression',fails)
  call shape_erm(x,kvals,er,'fh')
  call check(all(er%shape>0.0_dp).and.all(er%rho<0.0_dp),'FH exponential regression',fails)

  sig2=0.0_dp;sig2(1,1)=1.0_dp;sig2(2,2)=1.0_dp;up2=0.0_dp
  v=mvn_upper_prob_qmc(up2,sig2,4096)
  call check(abs(v-0.25_dp)<0.02_dp,'native bivariate normal probability',fails)
  z2=[1.0_dp,1.0_dp]
  call check(abs(expme_logistic(z2,0.5_dp)-sqrt(2.0_dp))<1.0e-14_dp,'logistic exponent measure',fails)
  lam2=reshape([0.0_dp,0.5_dp,0.5_dp,0.0_dp],[2,2])
  expect=2.0_dp*normal_cdf(sqrt(0.5_dp))
  call check(abs(expme_br(z2,lam2)-expect)<1.0e-12_dp,'Brown-Resnick bivariate exponent measure',fails)

  sig3=reshape([1.0_dp,0.3_dp,0.2_dp,0.3_dp,1.0_dp,0.25_dp,0.2_dp,0.25_dp,1.0_dp],[3,3])
  z3=[1.0_dp,1.2_dp,0.8_dp]
  v=expme_br_wt(z3,sig3,4096)
  call check(v>0.0_dp.and.v<sum(1.0_dp/z3)+0.2_dp,'BR-WT exponent measure finite',fails)
  v=expme_xstud(z3,sig3,4.0_dp,4096)
  call check(v>0.0_dp.and.v<sum(1.0_dp/z3)+0.2_dp,'extremal-t exponent measure finite',fails)

  qmat=0.0_dp;qmat(1,1)=2.0_dp;qmat(2,2)=2.0_dp;qmat(3,3)=2.0_dp
  qmat(1,2)=-0.5_dp;qmat(2,1)=-0.5_dp;qmat(1,3)=-0.4_dp;qmat(3,1)=-0.4_dp
  qmat(2,3)=-0.3_dp;qmat(3,2)=-0.3_dp;lvec=[0.1_dp,-0.2_dp,0.05_dp]
  v=expme_hr(z3,qmat,lvec,4096)
  call check(v>0.0_dp,'Huesler-Reiss exponent measure finite',fails)

  allocate(umat(n,2))
  do i=1,n
    umat(i,1)=(real(i,dp)-0.25_dp)/real(n+1,dp)
    umat(i,2)=min(0.999_dp,max(0.001_dp,umat(i,1)**0.92_dp))
  end do
  call kjtail_uniform(umat,[0.90_dp,0.95_dp],kj)
  call check(all(kj%n_tail>=3).and.all(kj%eta>0.0_dp).and.all(kj%eta<=1.0_dp), &
    'Krupskii-Joe tail dependence fit',fails)
  call check(all(kj%k1>0.0_dp),'Krupskii-Joe positive k1',fails)

  do i=1,n
    x(i)=qgp((real(i,dp)-0.5_dp)/real(n,dp),scale=1.2_dp,shape=0.10_dp)
  end do
  call gpd_profile(x,[-0.05_dp,0.10_dp,0.25_dp],'shape',pr)
  call check(all(pr%loglik>-huge(1.0_dp)/2.0_dp),'GPD profile likelihood finite',fails)
  call check(pr%loglik(2)>=pr%loglik(1)-2.0_dp.and.pr%loglik(2)>=pr%loglik(3)-2.0_dp, &
    'GPD profile peaks near generating shape',fails)

  allocate(xgev(n))
  do i=1,n
    xgev(i)=qgev((real(i,dp)-0.5_dp)/real(n,dp),loc=0.3_dp,scale=1.1_dp,shape=0.08_dp)
  end do
  call gev_profile(xgev,[-0.05_dp,0.08_dp,0.22_dp],'shape',pr)
  call check(all(pr%loglik>-huge(1.0_dp)/2.0_dp),'GEV profile likelihood finite',fails)
  call check(pr%loglik(2)>=pr%loglik(1)-2.0_dp.and.pr%loglik(2)>=pr%loglik(3)-2.0_dp, &
    'GEV profile peaks near generating shape',fails)

  if(fails==0)then
    print '(a)','test_parity_v02: PASS'
  else
    print '(a,i0)','test_parity_v02: FAIL ',fails
    error stop 1
  end if
contains
  subroutine check(ok,name,nfail)
    logical,intent(in)::ok
    character(len=*),intent(in)::name
    integer,intent(inout)::nfail
    if(.not.ok)then;print '(a,a)','FAIL: ',trim(name);nfail=nfail+1;end if
  end subroutine check
end program test_parity_v02
