program test_mev
  use mev_kinds, only: dp
  use mev_distributions
  use mev_univariate
  use mev_tailindex
  use mev_spatial
  use mev_sampling
  use mev_egp
  use mev_diagnostics
  use mev_emplik
  use mev_reparam
  implicit none
  integer :: fails, i, info
  real(dp) :: p, q, ll1, ll2, scale, shape
  real(dp) :: sampg(2500), sampgev(1200), dir(300,3), spec(100,3)
  real(dp) :: mx(1200,2), medx
  real(dp) :: alpha(3), mu3(3), sig3(3,3), mn(200,3)
  real(dp) :: locs(3,2), dm(3,3), lam(3,3), cov2(2,2)
  integer :: sub(2)
  real(dp) :: a(4), lm(4), thr(3), me(3), se(3)
  integer :: ne(3)
  real(dp) :: udat(6,2), qlev(1), eta(1), chi(1), cb(1)
  real(dp) :: z(4,1), target(1), ew(4), ang(3), wts(3), ss(2), pk(2)
  type(emplik_result) :: er
  real(dp) :: rls(50,3), sc3(3), inf3(3,3), pars3(3)
  type(mev_fit_result) :: gf, vf
  type(egp_fit_result) :: ef

  fails=0
  call seed_fixed()

  do i=1,5
    p=real(i,dp)/6.0_dp
    q=qgev(p,0.3_dp,1.7_dp,0.2_dp)
    call check(abs(pgev(q,0.3_dp,1.7_dp,0.2_dp)-p)<2.0e-10_dp,'GEV inversion',fails)
    q=qgp(p,scale=1.4_dp,shape=-0.15_dp)
    call check(abs(pgp(q,scale=1.4_dp,shape=-0.15_dp)-p)<2.0e-10_dp,'GPD inversion',fails)
  end do
  call check(abs(pgev(0.0_dp,0.0_dp,1.0_dp,0.0_dp)-exp(-1.0_dp))<1e-13_dp,'Gumbel identity',fails)

  do i=1,7
    p=0.73_dp
    q=qegp(p,1.3_dp,0.15_dp,1.4_dp,egp_name(i))
    call check(abs(pegp(q,1.3_dp,0.15_dp,1.4_dp,egp_name(i))-p)<5e-8_dp,'EGP inversion',fails)
  end do

  call rgp(size(sampg),sampg,scale=1.5_dp,shape=0.12_dp)
  call gpd_fit(sampg,gf)
  call check(gf%convergence==0,'GPD fit convergence',fails)
  call check(abs(gf%estimate(1)-1.5_dp)<0.25_dp,'GPD scale recovery',fails)
  call check(abs(gf%estimate(2)-0.12_dp)<0.12_dp,'GPD shape recovery',fails)

  call rgev(size(sampgev),sampgev,loc=0.4_dp,scale=1.1_dp,shape=-0.08_dp)
  call gev_fit(sampgev,vf)
  call check(vf%convergence==0,'GEV fit convergence',fails)
  call check(abs(vf%estimate(1)-0.4_dp)<0.18_dp,'GEV loc recovery',fails)

  call pwm([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp],a)
  call lmoments([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp],lm)
  call check(abs(a(1)-3.0_dp)<1e-14_dp,'PWM mean',fails)
  call check(abs(lm(1)-3.0_dp)<1e-14_dp,'L moment first',fails)
  call gpd_lmom(sampg,scale,shape)
  call check(scale>0.0_dp,'GPD Lmoment scale positive',fails)

  thr=[0.0_dp,1.0_dp,2.0_dp]
  call mrl_profile([0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp],thr,me,se,ne)
  call check(ne(2)==3,'MRL counts',fails)
  call check(abs(me(2)-2.0_dp)<1e-14_dp,'MRL mean',fails)

  locs=reshape([0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,2.0_dp],[3,2])
  call distg(locs,1.0_dp,0.0_dp,dm,info)
  call check(info==0.and.maxval(abs(dm-transpose(dm)))<1e-14_dp,'distg symmetry',fails)
  call check(abs(powerexp_cor(1.0_dp,1.0_dp,1.0_dp)-exp(-1.0_dp))<1e-14_dp,'powerexp',fails)
  lam=reshape([0.0_dp,0.25_dp,1.0_dp,0.25_dp,0.0_dp,0.36_dp,1.0_dp,0.36_dp,0.0_dp],[3,3])
  sub=[2,3];call lambda2cov(lam,1,sub,sub,cov2,info)
  call check(info==0.and.abs(cov2(1,1)-1.0_dp)<1e-14_dp,'Lambda2cov diagonal',fails)

  alpha=[0.5_dp,1.5_dp,2.0_dp];call rdir(300,alpha,dir,.true.,info)
  call check(info==0.and.maxval(abs(sum(dir,dim=2)-1.0_dp))<1e-12_dp,'Dirichlet simplex',fails)
  call rlogspec(100,3,2.0_dp,spec,info)
  call check(info==0.and.maxval(abs(sum(spec,dim=2)-1.0_dp))<1e-12_dp,'log spectral simplex',fails)
  call rneglogspec(100,3,1.5_dp,spec,info)
  call check(info==0.and.minval(spec)>=0.0_dp,'neglog spectral support',fails)
  call rmev(1200,'log',[2.0_dp],mx,info=info)
  call check(info==0.and.minval(mx)>0.0_dp,'rmev positive max-stable sample',fails)
  medx=median_column(mx(:,1))
  call check(abs(medx-1.0_dp/log(2.0_dp))<0.25_dp,'rmev unit Frechet margin',fails)

  mu3=0.0_dp;sig3=0.0_dp;sig3(1,1)=1.0_dp;sig3(2,2)=2.0_dp;sig3(3,3)=0.5_dp
  call rmnorm(200,mu3,sig3,mn,info)
  call check(info==0,'MV normal sample',fails)

  udat=reshape([0.90_dp,0.92_dp,0.70_dp,0.95_dp,0.85_dp,0.99_dp, &
                0.91_dp,0.93_dp,0.65_dp,0.96_dp,0.80_dp,0.98_dp],[6,2])
  qlev=[0.8_dp];call taildep_empirical(udat,qlev,eta,chi,cb,.true.)
  call check(chi(1)>=0.0_dp.and.eta(1)>=0.0_dp,'tail dependence empirical',fails)

  z(:,1)=[-1.5_dp,-0.5_dp,0.5_dp,1.5_dp];target=[0.0_dp]
  call euclidean_weights(z,target,ew,info)
  call check(info==0.and.abs(sum(ew)-1.0_dp)<1e-12_dp,'Euclidean weights',fails)
  call empirical_likelihood(z,target,er)
  call check(er%converged.and.abs(sum(er%weights)-1.0_dp)<1e-10_dp,'empirical likelihood',fails)
  ang=[0.2_dp,0.5_dp,0.8_dp];wts=[1.0_dp/3.0_dp,1.0_dp/3.0_dp,1.0_dp/3.0_dp];ss=[0.25_dp,0.5_dp]
  call pickands_empirical(ss,ang,wts,pk)
  call check(all(pk>=0.5_dp),'Pickands empirical',fails)

  call rrlarg(50,3,0.0_dp,1.0_dp,0.1_dp,rls)
  call check(all(rls(:,1)>=rls(:,2)).and.all(rls(:,2)>=rls(:,3)),'r-largest ordering',fails)
  pars3=[0.0_dp,1.0_dp,0.1_dp];ll1=rlarg_ll(pars3,rls)
  call rlarg_score(pars3,rls,sc3);call rlarg_infomat(pars3,rls,inf3)
  call check(ll1>-huge(1.0_dp)/2.0_dp.and.all(abs(sc3)<huge(1.0_dp)),'r-largest likelihood',fails)

  ll1=gpdr_ll([2.0_dp,0.2_dp],[0.2_dp,0.5_dp,1.0_dp],10.0_dp)
  scale=2.0_dp*0.2_dp/(10.0_dp**0.2_dp-1.0_dp)
  ll2=gpd_ll([scale,0.2_dp],[0.2_dp,0.5_dp,1.0_dp])
  call check(abs(ll1-ll2)<1e-12_dp,'GP return parameterization',fails)
  call check(pp_ll([0.0_dp,1.0_dp,0.0_dp],[1.0_dp,2.0_dp],0.5_dp,1.0_dp)<0.0_dp,'PP likelihood',fails)

  call regp(500,sampg(1:500),1.2_dp,0.1_dp,1.4_dp,'pt-power')
  call egp_fit(sampg(1:500),0.0_dp,'pt-power',ef,info=info)
  call check(info==0.and.ef%scale>0.0_dp.and.ef%kappa>0.0_dp,'EGP fit smoke',fails)

  if(fails==0)then
    print '(a)', 'test_mev: PASS'
  else
    print '(a,i0)', 'test_mev: FAIL ',fails
    error stop 1
  end if
contains
  subroutine check(ok,name,nfail)
    logical,intent(in)::ok
    character(len=*),intent(in)::name
    integer,intent(inout)::nfail
    if(.not.ok)then;print '(a,a)','FAIL: ',trim(name);nfail=nfail+1;end if
  end subroutine check

  subroutine seed_fixed()
    integer,allocatable::seed(:);integer::n,j
    call random_seed(size=n);allocate(seed(n));do j=1,n;seed(j)=137+37*j;end do;call random_seed(put=seed)
  end subroutine seed_fixed


  function median_column(x) result(m)
    real(dp),intent(in)::x(:)
    real(dp)::m,y(size(x)),t
    integer::ii,jj
    y=x
    do ii=2,size(y)
      t=y(ii);jj=ii-1
      do while(jj>=1)
        if(y(jj)<=t)exit
        y(jj+1)=y(jj);jj=jj-1
      end do
      y(jj+1)=t
    end do
    if(mod(size(y),2)==1)then
      m=y((size(y)+1)/2)
    else
      m=0.5_dp*(y(size(y)/2)+y(size(y)/2+1))
    end if
  end function median_column

  function egp_name(i) result(name)
    integer,intent(in)::i
    character(len=16)::name
    select case(i)
    case(1);name='pt-beta'
    case(2);name='pt-gamma'
    case(3);name='pt-power'
    case(4);name='gj-tnorm'
    case(5);name='gj-beta'
    case(6);name='exptilt'
    case default;name='logist'
    end select
  end function egp_name
end program test_mev
