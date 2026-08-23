module mev_profile
  use mev_kinds, only: dp
  use mev_math, only: pattern_minimize, variance_real, median_real
  use mev_univariate, only: gpd_ll, gev_ll
  use mev_reparam, only: gpdr_ll, gpde_ll, gpdn_ll, gevr_ll, gevn_ll
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: profile_result, gpd_profile, gev_profile

  type :: profile_result
    real(dp), allocatable :: psi(:)
    real(dp), allocatable :: loglik(:)
    real(dp), allocatable :: nuisance(:,:)
    integer, allocatable :: convergence(:)
  end type profile_result

contains

  subroutine gpd_profile(dat,psi,param,result,m,nblock)
    real(dp), intent(in) :: dat(:),psi(:)
    character(len=*), intent(in) :: param
    type(profile_result), intent(out) :: result
    real(dp), intent(in), optional :: m,nblock
    real(dp) :: x(1),fval,mm,nb,sc0,sh0,current_psi
    integer :: i,info
    character(len=16) :: pa
    pa=trim(adjustl(param)); mm=100.0_dp; if(present(m)) mm=m
    nb=100.0_dp; if(present(nblock)) nb=nblock
    allocate(result%psi(size(psi)),result%loglik(size(psi)),result%nuisance(size(psi),1),result%convergence(size(psi)))
    result%psi=psi;result%loglik=-huge(1.0_dp);result%nuisance=0.0_dp;result%convergence=1
    sc0=max(sqrt(variance_real(dat)),max(sum(dat)/max(1.0_dp,real(size(dat),dp)),1.0e-3_dp));sh0=0.1_dp
    do i=1,size(psi)
      current_psi=psi(i)
      select case(pa)
      case('shape')
        x(1)=log(sc0); call pattern_minimize(gpd_obj,x,fval,info,1500,1.0e-8_dp,0.2_dp)
        sc0=exp(x(1)); result%nuisance(i,1)=sc0
      case default
        x(1)=sh0; call pattern_minimize(gpd_obj,x,fval,info,1500,1.0e-8_dp,0.15_dp)
        sh0=x(1); result%nuisance(i,1)=sh0
      end select
      result%loglik(i)=-fval;result%convergence(i)=info
    end do
  contains
    real(dp) function gpd_obj(y) result(v)
      real(dp),intent(in)::y(:)
      real(dp)::ll,sc,sh
      select case(pa)
      case('scale')
        sc=current_psi;sh=y(1);ll=gpd_ll([sc,sh],dat)
      case('shape')
        sc=exp(min(50.0_dp,max(-50.0_dp,y(1))));sh=current_psi;ll=gpd_ll([sc,sh],dat)
      case('quant','VaR')
        ll=gpdr_ll([current_psi,y(1)],dat,mm)
      case('ES','es')
        ll=gpde_ll([current_psi,y(1)],dat,mm)
      case('Nmean','nmean')
        ll=gpdn_ll([current_psi,y(1)],dat,nb)
      case default
        ll=-huge(1.0_dp)
      end select
      if(.not.ieee_is_finite(ll)) then;v=1.0e100_dp;else;v=-ll;end if
    end function gpd_obj
  end subroutine gpd_profile

  subroutine gev_profile(dat,psi,param,result,p,nblock,q)
    real(dp), intent(in) :: dat(:),psi(:)
    character(len=*), intent(in) :: param
    type(profile_result), intent(out) :: result
    real(dp), intent(in), optional :: p,nblock,q
    real(dp) :: x(2),fval,pp,nb,qq,mu0,sc0,sh0,current_psi
    integer :: i,info
    character(len=16) :: pa
    pa=trim(adjustl(param));pp=0.01_dp;if(present(p))pp=p;nb=100.0_dp;if(present(nblock))nb=nblock
    qq=0.5_dp;if(present(q))qq=q
    allocate(result%psi(size(psi)),result%loglik(size(psi)),result%nuisance(size(psi),2),result%convergence(size(psi)))
    result%psi=psi;result%loglik=-huge(1.0_dp);result%nuisance=0.0_dp;result%convergence=1
    mu0=median_real(dat);sc0=max(sqrt(variance_real(dat))*0.8_dp,1.0e-3_dp);sh0=0.1_dp
    do i=1,size(psi)
      current_psi=psi(i)
      select case(pa)
      case('loc')
        x=[log(sc0),sh0]
      case('scale')
        x=[mu0,sh0]
      case('shape')
        x=[mu0,log(sc0)]
      case('quant','VaR')
        x=[log(sc0),sh0]
      case('Nmean','nmean','Nquant','nquant')
        x=[mu0,sh0]
      case default
        cycle
      end select
      call pattern_minimize(gev_obj,x,fval,info,2200,1.0e-8_dp,0.18_dp)
      select case(pa)
      case('loc','quant','VaR')
        sc0=exp(x(1));sh0=x(2);result%nuisance(i,:)=[sc0,sh0]
      case('scale')
        mu0=x(1);sh0=x(2);result%nuisance(i,:)=[mu0,sh0]
      case('shape')
        mu0=x(1);sc0=exp(x(2));result%nuisance(i,:)=[mu0,sc0]
      case default
        mu0=x(1);sh0=x(2);result%nuisance(i,:)=[mu0,sh0]
      end select
      result%loglik(i)=-fval;result%convergence(i)=info
    end do
  contains
    real(dp) function gev_obj(y) result(v)
      real(dp),intent(in)::y(:)
      real(dp)::ll,mu,sc,sh
      select case(pa)
      case('loc')
        mu=current_psi;sc=exp(min(50.0_dp,max(-50.0_dp,y(1))));sh=y(2);ll=gev_ll([mu,sc,sh],dat)
      case('scale')
        mu=y(1);sc=current_psi;sh=y(2);ll=gev_ll([mu,sc,sh],dat)
      case('shape')
        mu=y(1);sc=exp(min(50.0_dp,max(-50.0_dp,y(2))));sh=current_psi;ll=gev_ll([mu,sc,sh],dat)
      case('quant','VaR')
        sc=exp(min(50.0_dp,max(-50.0_dp,y(1))));sh=y(2);ll=gevr_ll([current_psi,sc,sh],dat,pp)
      case('Nmean','nmean')
        ll=gevn_ll([y(1),current_psi,y(2)],dat,nb,qq,'mean')
      case('Nquant','nquant')
        ll=gevn_ll([y(1),current_psi,y(2)],dat,nb,qq,'quantile')
      case default
        ll=-huge(1.0_dp)
      end select
      if(.not.ieee_is_finite(ll)) then;v=1.0e100_dp;else;v=-ll;end if
    end function gev_obj
  end subroutine gev_profile

end module mev_profile
