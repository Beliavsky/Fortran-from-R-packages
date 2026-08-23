module mev_univariate
  use mev_kinds, only: dp, pi
  use mev_math, only: digamma_mev, finite_diff_hessian, pattern_minimize, &
      mean_real, variance_real, sort_ascending
  use mev_distributions, only: qgev
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
  implicit none
  private
  public :: mev_fit_result, gpd_ll, gpd_score, gpd_infomat, gev_ll, gev_score, gev_infomat
  public :: gpd_fit, gev_fit, gev_retlev, gpd_n_mean, gpd_n_quant, gev_n_mean, gev_n_quant
  public :: gev_nyr, gpd_to_pareto, jac_gpd_pareto, ordexp_to_gev

  type :: mev_fit_result
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: vcov(:,:)
    real(dp) :: loglik = -huge(1.0_dp)
    integer :: convergence = 1
    integer :: nobs = 0
  end type mev_fit_result

contains

  pure real(dp) function gpd_ll(par,dat,tol) result(ll)
    real(dp),intent(in)::par(2),dat(:)
    real(dp),intent(in),optional::tol
    real(dp)::sigma,xi,toluse,t(size(dat))
    integer::n
    toluse=1.0e-5_dp;if(present(tol))toluse=tol
    sigma=par(1);xi=par(2);n=size(dat)
    if(sigma<=0.0_dp.or.any(dat<0.0_dp))then
      ll=-ieee_value(0.0_dp,ieee_positive_inf);return
    end if
    if(abs(xi)>toluse)then
      if(xi>-1.0_dp)then
        t=1.0_dp+xi*dat/sigma
        if(any(t<=0.0_dp))then
          ll=-ieee_value(0.0_dp,ieee_positive_inf)
        else
          ll=-real(n,dp)*log(sigma)-(1.0_dp+1.0_dp/xi)*sum(log(t))
        end if
      else if(abs(xi+1.0_dp)<1.0e-12_dp)then
        if(maxval(dat)<=sigma)then;ll=-real(n,dp)*log(sigma);else;ll=-ieee_value(0.0_dp,ieee_positive_inf);end if
      else
        ll=-ieee_value(0.0_dp,ieee_positive_inf)
      end if
    else
      ll=-real(n,dp)*log(sigma)-sum(dat)/sigma
    end if
  end function gpd_ll

  pure subroutine gpd_score(par,dat,score)
    real(dp),intent(in)::par(2),dat(:)
    real(dp),intent(out)::score(2)
    real(dp)::sigma,xi,t(size(dat))
    sigma=par(1);xi=par(2)
    if(abs(xi)>1.0e-10_dp)then
      t=1.0_dp+xi*dat/sigma
      if(any(t<=0.0_dp).or.sigma<=0.0_dp)then
        score=ieee_value(0.0_dp,ieee_quiet_nan);return
      end if
      score(1)=sum(dat*(xi+1.0_dp)/(sigma*sigma*t)-1.0_dp/sigma)
      score(2)=sum(-dat*(1.0_dp/xi+1.0_dp)/(sigma*t)+log(t)/(xi*xi))
    else
      score(1)=sum((dat-sigma)/(sigma*sigma))
      score(2)=sum(0.5_dp*(dat-2.0_dp*sigma)*dat/(sigma*sigma))
    end if
  end subroutine gpd_score

  subroutine gpd_infomat(par,dat,info,expected,nobs)
    real(dp),intent(in)::par(2),dat(:)
    real(dp),intent(out)::info(2,2)
    logical,intent(in),optional::expected
    integer,intent(in),optional::nobs
    real(dp)::sigma,xi,k11,k12,k22,h(2,2)
    integer::nn
    logical::ex
    ex=.false.;if(present(expected))ex=expected
    sigma=par(1);xi=par(2)
    if(ex)then
      nn=size(dat);if(present(nobs))nn=nobs
      if(xi<=-0.5_dp.or.sigma<=0.0_dp)then
        info=ieee_value(0.0_dp,ieee_quiet_nan);return
      end if
      k22=-2.0_dp/((1.0_dp+xi)*(1.0_dp+2.0_dp*xi))
      k11=-1.0_dp/(sigma*sigma*(1.0_dp+2.0_dp*xi))
      k12=-1.0_dp/(sigma*(1.0_dp+xi)*(1.0_dp+2.0_dp*xi))
      info=-real(nn,dp)*reshape([k11,k12,k12,k22],[2,2])
    else
      call finite_diff_hessian(nll_local,par,h,1.0e-4_dp)
      info=h
    end if
  contains
    function nll_local(x) result(v)
      real(dp),intent(in)::x(:)
      real(dp)::v
      v=-gpd_ll(x(1:2),dat)
      if(.not.ieee_is_finite(v))v=huge(1.0_dp)/100.0_dp
    end function
  end subroutine gpd_infomat

  pure real(dp) function gev_ll(par,dat) result(ll)
    real(dp),intent(in)::par(3),dat(:)
    real(dp)::mu,sigma,xi,t(size(dat)),z(size(dat))
    mu=par(1);sigma=par(2);xi=par(3)
    if(sigma<=0.0_dp)then;ll=-ieee_value(0.0_dp,ieee_positive_inf);return;end if
    z=(dat-mu)/sigma
    if(abs(xi)>1.0e-7_dp)then
      t=1.0_dp+xi*z
      if(any(t<=0.0_dp))then
        ll=-ieee_value(0.0_dp,ieee_positive_inf)
      else
        ll=sum(-log(sigma)-(1.0_dp/xi+1.0_dp)*log(t)-t**(-1.0_dp/xi))
      end if
    else
      ll=sum(-log(sigma)-z-exp(-z))
    end if
  end function gev_ll

  pure subroutine gev_score(par,dat,score)
    real(dp),intent(in)::par(3),dat(:)
    real(dp),intent(out)::score(3)
    real(dp)::mu,sigma,xi,z(size(dat)),t(size(dat)),lt(size(dat)),a(size(dat))
    mu=par(1);sigma=par(2);xi=par(3)
    z=(dat-mu)/sigma
    if(sigma<=0.0_dp)then;score=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    if(abs(xi)>1.0e-10_dp)then
      t=1.0_dp+xi*z
      if(any(t<=0.0_dp))then;score=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
      lt=log(t);a=t**(-1.0_dp/xi)
      score(1)=sum((1.0_dp+xi-a)/(sigma*t))
      score(2)=sum(-1.0_dp/sigma + z*(1.0_dp+xi-a)/(sigma*t))
      score(3)=sum(lt/(xi*xi)*(1.0_dp-a) + z/(xi*t)*(a-1.0_dp-xi))
    else
      score(1)=sum((1.0_dp-exp(-z))/sigma)
      score(2)=sum((-1.0_dp+z*(1.0_dp-exp(-z)))/sigma)
      score(3)=sum(0.5_dp*z*(exp(-z)*(-z)+z-2.0_dp))
    end if
  end subroutine gev_score

  subroutine gev_infomat(par,dat,info,expected,nobs)
    real(dp),intent(in)::par(3),dat(:)
    real(dp),intent(out)::info(3,3)
    logical,intent(in),optional::expected
    integer,intent(in),optional::nobs
    real(dp)::sigma,xi,pv,qv,g2,h(3,3),psi1
    integer::nn
    logical::ex
    ex=.false.;if(present(expected))ex=expected
    sigma=par(2);xi=par(3)
    if(.not.ex)then
      call finite_diff_hessian(nll_local,par,h,1.0e-4_dp)
      info=h
      return
    end if
    nn=size(dat);if(present(nobs))nn=nobs
    if(sigma<=0.0_dp.or.xi< -0.5_dp)then;info=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    if(abs(xi)<1.0e-3_dp)then
      info=real(nn,dp)*reshape([ &
        1.0_dp/sigma**2, -0.422784335098467_dp/sigma**2, 0.41184033042644_dp/sigma, &
        -0.422784335098467_dp/sigma**2, 1.82368066085288_dp/sigma**2, 0.332484907160274_dp/sigma, &
        0.41184033042644_dp/sigma, 0.332484907160274_dp/sigma, 2.42360605517703_dp],[3,3])
      return
    end if
    pv=(1.0_dp+xi)**2*gamma(1.0_dp+2.0_dp*xi)
    qv=gamma(2.0_dp+xi)*(digamma_mev(1.0_dp+xi)+(1.0_dp+xi)/xi)
    psi1=digamma_mev(1.0_dp)
    g2=gamma(2.0_dp+xi)
    info(1,1)=pv/sigma**2
    info(1,2)=-(pv-g2)/(sigma**2*xi)
    info(1,3)=(pv/xi-qv)/(sigma*xi)
    info(2,1)=info(1,2)
    info(2,2)=(1.0_dp-2.0_dp*g2+pv)/(sigma**2*xi**2)
    info(2,3)=-(1.0_dp+psi1+(1.0_dp-g2)/xi-qv+pv/xi)/(sigma*xi**2)
    info(3,1)=info(1,3);info(3,2)=info(2,3)
    info(3,3)=(pi*pi/6.0_dp+(1.0_dp+psi1+1.0_dp/xi)**2-2.0_dp*qv/xi+pv/xi**2)/xi**2
    info=real(nn,dp)*info
  contains
    function nll_local(x) result(v)
      real(dp),intent(in)::x(:)
      real(dp)::v
      v=-gev_ll(x(1:3),dat)
      if(.not.ieee_is_finite(v))v=huge(1.0_dp)/100.0_dp
    end function
  end subroutine gev_infomat

  subroutine gpd_fit(dat,result)
    real(dp),intent(in)::dat(:)
    type(mev_fit_result),intent(out)::result
    real(dp)::theta(2),fval,par(2),info(2,2),inv(2,2),m,v,xi0,sig0
    integer::ier
    if(size(dat)<2.or.any(dat<0.0_dp))then
      allocate(result%estimate(2),result%vcov(2,2));result%estimate=0;result%vcov=0;result%convergence=2;return
    end if
    m=mean_real(dat);v=variance_real(dat)
    xi0=0.5_dp*(1.0_dp-m*m/max(v,tiny(1.0_dp)))
    xi0=max(-0.4_dp,min(0.5_dp,xi0))
    sig0=max(tiny(1.0_dp),m*(1.0_dp-xi0))
    theta=[log(sig0),xi0]
    call pattern_minimize(obj,theta,fval,ier,maxiter=1500,tol=1.0e-8_dp,initial_step=0.2_dp)
    par=[exp(theta(1)),theta(2)]
    allocate(result%estimate(2),result%vcov(2,2));result%estimate=par
    result%loglik=-fval;result%convergence=ier;result%nobs=size(dat)
    call gpd_infomat(par,dat,info,.false.)
    call inverse2(info,inv,ier);if(ier==0)then;result%vcov=inv;else;result%vcov=0.0_dp;end if
  contains
    function obj(th) result(vv)
      real(dp),intent(in)::th(:)
      real(dp)::vv,pp(2)
      pp=[exp(th(1)),th(2)]
      if(pp(2)<=-1.0_dp)then;vv=huge(1.0_dp)/100.0_dp;return;end if
      vv=-gpd_ll(pp,dat)
      if(.not.ieee_is_finite(vv))vv=huge(1.0_dp)/100.0_dp
    end function
  end subroutine gpd_fit

  subroutine gev_fit(dat,result)
    real(dp),intent(in)::dat(:)
    type(mev_fit_result),intent(out)::result
    real(dp)::theta(3),fval,par(3),info(3,3),inv(3,3),m,s
    integer::ier
    if(size(dat)<3)then
      allocate(result%estimate(3),result%vcov(3,3));result%estimate=0;result%vcov=0;result%convergence=2;return
    end if
    m=mean_real(dat);s=sqrt(max(variance_real(dat),tiny(1.0_dp)))
    theta=[m-0.5772156649015329_dp*s*sqrt(6.0_dp)/pi,log(max(s*sqrt(6.0_dp)/pi,tiny(1.0_dp))),0.0_dp]
    call pattern_minimize(obj,theta,fval,ier,maxiter=2500,tol=5.0e-8_dp,initial_step=0.15_dp)
    par=[theta(1),exp(theta(2)),theta(3)]
    allocate(result%estimate(3),result%vcov(3,3));result%estimate=par
    result%loglik=-fval;result%convergence=ier;result%nobs=size(dat)
    call gev_infomat(par,dat,info,.false.)
    call inverse3(info,inv,ier);if(ier==0)then;result%vcov=inv;else;result%vcov=0.0_dp;end if
  contains
    function obj(th) result(vv)
      real(dp),intent(in)::th(:)
      real(dp)::vv,pp(3)
      pp=[th(1),exp(th(2)),th(3)]
      if(pp(3)<=-0.99_dp)then;vv=huge(1.0_dp)/100.0_dp;return;end if
      vv=-gev_ll(pp,dat)
      if(.not.ieee_is_finite(vv))vv=huge(1.0_dp)/100.0_dp
    end function
  end subroutine gev_fit

  pure real(dp) function gev_retlev(par,p) result(v)
    real(dp),intent(in)::par(3),p
    real(dp)::yp
    yp=-log(1.0_dp-p)
    if(abs(par(3))<1.0e-10_dp)then
      v=par(1)-par(2)*log(yp)
    else
      v=par(1)-par(2)/par(3)*(1.0_dp-yp**(-par(3)))
    end if
  end function gev_retlev

  pure real(dp) function gpd_n_mean(par,nblock) result(v)
    real(dp),intent(in)::par(2)
    integer,intent(in)::nblock
    real(dp)::xi,sig,nn
    sig=par(1);xi=par(2);nn=real(nblock,dp)
    if(abs(xi)>1.0e-10_dp)then
      if(xi>=1.0_dp)then;v=ieee_value(0.0_dp,ieee_positive_inf);else
        v=(exp(log_gamma(nn+1.0_dp)+log_gamma(1.0_dp-xi)-log_gamma(nn+1.0_dp-xi))-1.0_dp)*sig/xi
      end if
    else
      v=(digamma_mev(nn+1.0_dp)-digamma_mev(1.0_dp))*sig
    end if
  end function

  pure real(dp) function gpd_n_quant(par,q,nblock) result(v)
    real(dp),intent(in)::par(2),q
    integer,intent(in)::nblock
    if(abs(par(2))>1.0e-10_dp)then
      v=par(1)/par(2)*((1.0_dp-q**(1.0_dp/real(nblock,dp)))**(-par(2))-1.0_dp)
    else
      v=-par(1)*log(1.0_dp-q**(1.0_dp/real(nblock,dp)))
    end if
  end function

  pure real(dp) function gev_n_mean(par,nblock) result(v)
    real(dp),intent(in)::par(3)
    integer,intent(in)::nblock
    real(dp)::nn
    nn=real(nblock,dp)
    if(abs(par(3))>1.0e-10_dp)then
      if(par(3)>=1.0_dp)then;v=ieee_value(0.0_dp,ieee_positive_inf);else
        v=par(1)-par(2)/par(3)*(1.0_dp-nn**par(3)*gamma(1.0_dp-par(3)))
      end if
    else
      v=par(1)+par(2)*(log(nn)-digamma_mev(1.0_dp))
    end if
  end function

  pure real(dp) function gev_n_quant(par,nblock,q) result(v)
    real(dp),intent(in)::par(3),q
    integer,intent(in)::nblock
    real(dp)::nn
    nn=real(nblock,dp)
    if(abs(par(3))>1.0e-10_dp)then
      v=par(1)-par(2)/par(3)*(1.0_dp-(nn/log(1.0_dp/q))**par(3))
    else
      v=par(1)+par(2)*(log(nn)-log(log(1.0_dp/q)))
    end if
  end function

  subroutine gev_nyr(par,nobs,nblock,kind,p,estimate,variance)
    real(dp),intent(in)::par(3)
    integer,intent(in)::nobs,nblock
    character(len=*),intent(in)::kind
    real(dp),intent(in),optional::p
    real(dp),intent(out)::estimate,variance
    real(dp)::grad(3),info(3,3),inv(3,3),yp,nn,xi,sig,emu
    real(dp)::dummy(1)
    integer::ier
    dummy=0.0_dp;xi=par(3);sig=par(2);nn=real(nblock,dp);emu=-digamma_mev(1.0_dp)
    select case(trim(kind))
    case('retlev')
      if(.not.present(p))then;estimate=ieee_value(0.0_dp,ieee_quiet_nan);variance=estimate;return;end if
      yp=-log(1.0_dp-p)
      estimate=gev_retlev(par,p)
      if(abs(xi)<1.0e-10_dp)then
        grad=[1.0_dp,-log(yp),0.5_dp*sig*log(yp)**2]
      else
        grad=[1.0_dp,-(1.0_dp-yp**(-xi))/xi, &
          sig*(1.0_dp-yp**(-xi))/xi**2-sig/xi*yp**(-xi)*log(yp)]
      end if
    case('median')
      estimate=gev_n_quant(par,nblock,0.5_dp)
      if(abs(xi)<1.0e-10_dp)then
        grad=[1.0_dp,log(nn/log(2.0_dp)),0.5_dp*sig*log(nn/log(2.0_dp))**2]
      else
        grad=[1.0_dp,((nn/log(2.0_dp))**xi-1.0_dp)/xi, &
          sig*(nn/log(2.0_dp))**xi*log(nn/log(2.0_dp))/xi &
          -sig*((nn/log(2.0_dp))**xi-1.0_dp)/xi**2]
      end if
    case('mean')
      estimate=gev_n_mean(par,nblock)
      if(abs(xi)<1.0e-10_dp)then
        grad=[1.0_dp,log(nn)+emu,0.5_dp*sig*(emu**2+pi*pi/6.0_dp+2.0_dp*emu*log(nn)+log(nn)**2)]
      else
        grad=[1.0_dp,(nn**xi*gamma(1.0_dp-xi)-1.0_dp)/xi, &
          (nn**xi*log(nn)*gamma(1.0_dp-xi)-nn**xi*digamma_mev(1.0_dp-xi)*gamma(1.0_dp-xi))*sig/xi &
          -(nn**xi*gamma(1.0_dp-xi)-1.0_dp)*sig/xi**2]
      end if
    case default
      estimate=ieee_value(0.0_dp,ieee_quiet_nan);variance=estimate;return
    end select
    call gev_infomat(par,dummy,info,.true.,nobs)
    call inverse3(info,inv,ier)
    if(ier==0)then;variance=dot_product(grad,matmul(inv,grad));else;variance=ieee_value(0.0_dp,ieee_quiet_nan);end if
  end subroutine gev_nyr

  pure real(dp) function gpd_to_pareto(x,loc,scale,shape,lambdau) result(v)
    real(dp),intent(in)::x,loc,scale,shape,lambdau
    if(abs(shape)<1.0e-10_dp)then
      v=exp(max(x-loc,0.0_dp)/scale)/lambdau
    else
      v=max(0.0_dp,1.0_dp+shape/scale*max(x-loc,0.0_dp))**(1.0_dp/shape)/lambdau
    end if
  end function

  pure real(dp) function jac_gpd_pareto(dat,censored,loc,scale,shape,lambdau) result(ll)
    real(dp),intent(in)::dat(:),loc,scale,shape,lambdau
    logical,intent(in)::censored(:)
    integer::j,nunc
    ll=0.0_dp;nunc=0
    do j=1,size(dat)
      if(.not.censored(j))then
        nunc=nunc+1
        if(abs(shape)<1.0e-10_dp)then
          ll=ll+(dat(j)-loc)/scale
        else
          ll=ll+(1.0_dp/shape-1.0_dp)*log(1.0_dp+shape/scale*max(0.0_dp,dat(j)-loc))
        end if
      end if
    end do
    ll=ll-real(nunc,dp)*(log(scale)+log(lambdau))
  end function

  pure real(dp) function ordexp_to_gev(edat,par) result(v)
    real(dp),intent(in)::edat,par(3)
    if(abs(par(3))>1.0e-8_dp)then
      v=par(1)+par(2)*(edat**(-par(3))-1.0_dp)/par(3)
    else
      v=par(1)-par(2)*log(edat)
    end if
  end function

  subroutine inverse2(a,ainv,info)
    real(dp),intent(in)::a(2,2)
    real(dp),intent(out)::ainv(2,2)
    integer,intent(out)::info
    real(dp)::det
    det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
    if(abs(det)<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a))))then;info=1;ainv=0.0_dp;return;end if
    ainv=reshape([a(2,2),-a(2,1),-a(1,2),a(1,1)],[2,2])/det;info=0
  end subroutine

  subroutine inverse3(a,ainv,info)
    real(dp),intent(in)::a(3,3)
    real(dp),intent(out)::ainv(3,3)
    integer,intent(out)::info
    real(dp)::det
    det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
       -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
       +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    if(abs(det)<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a))))then;info=1;ainv=0.0_dp;return;end if
    ainv(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
    ainv(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
    ainv(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
    ainv(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
    ainv(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
    ainv(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
    ainv(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
    ainv(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
    ainv(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
    ainv=ainv/det;info=0
  end subroutine

end module mev_univariate
