! SPDX-License-Identifier: GPL-3.0-only
module mass_regression
  use rrcov_kinds, only : dp
  use rrcov_types, only : rrcov_success
  use rrcov_linalg, only : general_inverse, symmetric_eigen, identity_matrix
  use mass_types, only : regression_result, ridge_result, huber_result, mass_success, &
    mass_invalid_argument, mass_dimension_error, mass_no_convergence
  use mass_math, only : least_squares, weighted_least_squares, mad_mass, sample_sd, normal_pdf
  use rrcov_stats, only : normal_cdf
  implicit none
  private
  public :: linear_model_fit, lm_gls, lm_ridge, ridge_coefficients
  public :: psi_huber, psi_hampel, psi_bisquare, rlm_fit
  public :: standardized_residuals, studentized_residuals
  public :: boxcox_profile, logtrans_profile, dose_p
  public :: huber_location, huber_location_scale

contains

  subroutine linear_model_fit(x,y,result,weights)
    real(dp),intent(in)::x(:,:),y(:)
    type(regression_result),intent(out)::result
    real(dp),intent(in),optional::weights(:)
    real(dp),allocatable::inv(:,:),xtx(:,:),cov(:,:),w(:)
    real(dp)::rss
    integer::rank,status,i,n,p,st
    n=size(x,1); p=size(x,2)
    if(n/=size(y) .or. n<=p .or. p<1) then; result%status=mass_dimension_error; return; end if
    if(present(weights)) then
      call weighted_least_squares(x,y,weights,result%coefficients,result%residuals,cov,rank,status)
      w=weights
    else
      call least_squares(x,y,result%coefficients,result%residuals,rank,status)
      allocate(w(n)); w=1.0_dp
      xtx=matmul(transpose(x),x); inv=general_inverse(xtx,st); cov=inv
    end if
    result%rank=rank; result%df_residual=n-rank; result%fitted=y-result%residuals
    rss=sum(w*result%residuals**2)
    result%sigma=sqrt(rss/real(max(1,result%df_residual),dp))
    result%covariance=result%sigma**2*cov
    allocate(result%weights(n),result%leverages(n)); result%weights=w
    do i=1,n
      result%leverages(i)=w(i)*dot_product(x(i,:),matmul(cov,x(i,:)))
    end do
    result%log_likelihood=-0.5_dp*real(n,dp)*(log(2.0_dp*acos(-1.0_dp)*rss/real(n,dp))+1.0_dp)
    result%aic=-2.0_dp*result%log_likelihood+2.0_dp*real(rank+1,dp)
    result%status=status; result%method="ordinary least squares"
  end subroutine linear_model_fit

  subroutine lm_gls(x,y,wmat,result,inverse)
    real(dp),intent(in)::x(:,:),y(:),wmat(:,:)
    type(regression_result),intent(out)::result
    logical,intent(in),optional::inverse
    real(dp),allocatable::ev(:),vec(:,:),a(:,:),ainv(:,:),d(:,:),di(:,:),xt(:,:),yt(:),fit_t(:)
    logical::invflag
    integer::n,i,st
    n=size(y); invflag=.false.; if(present(inverse)) invflag=inverse
    if(size(x,1)/=n .or. size(wmat,1)/=n .or. size(wmat,2)/=n) then; result%status=mass_dimension_error; return; end if
    call symmetric_eigen(wmat,ev,vec,st)
    if(minval(ev)<=0.0_dp) then; result%status=mass_invalid_argument; return; end if
    allocate(d(n,n),di(n,n)); d=0.0_dp; di=0.0_dp
    do i=1,n
      if(invflag) then
        d(i,i)=ev(i)**(-0.5_dp); di(i,i)=ev(i)**0.5_dp
      else
        d(i,i)=ev(i)**0.5_dp; di(i,i)=ev(i)**(-0.5_dp)
      end if
    end do
    a=matmul(d,transpose(vec)); ainv=matmul(vec,di)
    xt=matmul(a,x); yt=matmul(a,y)
    call linear_model_fit(xt,yt,result)
    fit_t=result%fitted
    result%fitted=matmul(ainv,fit_t)
    result%residuals=y-result%fitted
    result%method="generalized least squares"
  end subroutine lm_gls

  subroutine lm_ridge(x,y,lambdas,result,intercept)
    real(dp),intent(in)::x(:,:),y(:),lambdas(:)
    type(ridge_result),intent(out)::result
    logical,intent(in),optional::intercept
    real(dp),allocatable::xc(:,:),yc(:),xm(:),sc(:),xtx(:,:),rhs(:),a(:,:),inv(:,:),fit(:)
    real(dp)::ym,rss,traceh
    logical::inter
    integer::n,p,k,j,st
    n=size(x,1); p=size(x,2); k=size(lambdas); inter=.true.; if(present(intercept)) inter=intercept
    if(n/=size(y) .or. n<=p .or. p<1 .or. k<1 .or. any(lambdas<0.0_dp)) then
      result%status=mass_invalid_argument; return
    end if
    allocate(xc(n,p),yc(n),xm(p),sc(p),result%coefficients(p,k),result%gcv(k),result%scales(p),result%lambdas(k))
    if(inter) then; xm=sum(x,dim=1)/real(n,dp); ym=sum(y)/real(n,dp); else; xm=0.0_dp; ym=0.0_dp; end if
    xc=x-spread(xm,1,n); yc=y-ym
    do j=1,p
      sc(j)=sqrt(sum(xc(:,j)**2)/real(n,dp)); if(sc(j)<=tiny(1.0_dp)) sc(j)=1.0_dp
      xc(:,j)=xc(:,j)/sc(j)
    end do
    xtx=matmul(transpose(xc),xc); rhs=matmul(transpose(xc),yc)
    do j=1,k
      a=xtx+lambdas(j)*identity_matrix(p); inv=general_inverse(a,st)
      result%coefficients(:,j)=matmul(inv,rhs)
      fit=matmul(xc,result%coefficients(:,j)); rss=sum((yc-fit)**2)
      traceh=sum([(dot_product(xtx(:,p),inv(p,:)),p=1,size(xtx,1))])
      result%gcv(j)=rss/max((real(n,dp)-traceh)**2,tiny(1.0_dp))
    end do
    result%scales=sc; result%lambdas=lambdas; result%y_mean=ym; result%status=mass_success
  end subroutine lm_ridge

  subroutine ridge_coefficients(result,x_mean,coefficients)
    type(ridge_result),intent(in)::result
    real(dp),intent(in)::x_mean(:)
    real(dp),allocatable,intent(out)::coefficients(:,:)
    integer::j,k,p
    p=size(result%coefficients,1); k=size(result%coefficients,2)
    allocate(coefficients(p+1,k))
    do j=1,k
      coefficients(2:,j)=result%coefficients(:,j)/result%scales
      coefficients(1,j)=result%y_mean-dot_product(coefficients(2:,j),x_mean)
    end do
  end subroutine ridge_coefficients

  pure elemental function psi_huber(u,k,deriv) result(value)
    real(dp),intent(in)::u
    real(dp),intent(in),optional::k
    integer,intent(in),optional::deriv
    real(dp)::value,kk
    integer::d
    kk=1.345_dp; if(present(k)) kk=k; d=0; if(present(deriv)) d=deriv
    if(d==0) then
      if(abs(u)<=kk) then; value=u; else; value=sign(kk,u); end if
    else
      value=merge(1.0_dp,0.0_dp,abs(u)<=kk)
    end if
  end function psi_huber

  pure elemental function psi_hampel(u,a,b,c,deriv) result(value)
    real(dp),intent(in)::u
    real(dp),intent(in),optional::a,b,c
    integer,intent(in),optional::deriv
    real(dp)::value,aa,bb,cc,z
    integer::d
    aa=2.0_dp;bb=4.0_dp;cc=8.0_dp
    if(present(a))aa=a;if(present(b))bb=b;if(present(c))cc=c;d=0;if(present(deriv))d=deriv;z=abs(u)
    if(d==0) then
      if(z<=aa) then; value=u
      else if(z<=bb) then; value=sign(aa,u)
      else if(z<=cc) then; value=sign(aa*(cc-z)/(cc-bb),u)
      else; value=0.0_dp; end if
    else
      if(z<=aa) then; value=1.0_dp
      else if(z<=bb) then; value=0.0_dp
      else if(z<=cc) then; value=-aa/(cc-bb)
      else; value=0.0_dp; end if
    end if
  end function psi_hampel

  pure elemental function psi_bisquare(u,c,deriv) result(value)
    real(dp),intent(in)::u
    real(dp),intent(in),optional::c
    integer,intent(in),optional::deriv
    real(dp)::value,cc,z,t
    integer::d
    cc=4.685_dp;if(present(c))cc=c;d=0;if(present(deriv))d=deriv;z=u/cc
    if(abs(z)>=1.0_dp) then; value=0.0_dp
    else
      t=1.0_dp-z*z
      if(d==0) then; value=u*t*t; else; value=t*(1.0_dp-5.0_dp*z*z); end if
    end if
  end function psi_bisquare

  subroutine rlm_fit(x,y,result,method,psi,scale_method,maxit,tolerance,case_weights)
    real(dp),intent(in)::x(:,:),y(:)
    type(regression_result),intent(out)::result
    character(len=*),intent(in),optional::method,psi,scale_method
    integer,intent(in),optional::maxit
    real(dp),intent(in),optional::tolerance,case_weights(:)
    character(len=16)::meth,psifun,scalem
    real(dp),allocatable::beta(:),oldbeta(:),res(:),w(:),cov(:,:),cw(:)
    real(dp)::s,tol,u,delta
    integer::rank,status,it,mit,i,n
    meth="M";if(present(method))meth=method;psifun="huber";if(present(psi))psifun=psi
    scalem="MAD";if(present(scale_method))scalem=scale_method;mit=50;if(present(maxit))mit=maxit
    tol=1.0e-4_dp;if(present(tolerance))tol=tolerance;n=size(y)
    if(size(x,1)/=n .or. n<=size(x,2)) then; result%status=mass_dimension_error; return; end if
    allocate(cw(n));cw=1.0_dp;if(present(case_weights)) then
      if(size(case_weights)/=n .or. any(case_weights<0.0_dp)) then; result%status=mass_invalid_argument; return; end if
      cw=case_weights
    end if
    call weighted_least_squares(x,y,cw,beta,res,cov,rank,status)
    s=mad_mass(res);if(s<=tiny(1.0_dp))s=max(sample_sd(res),sqrt(epsilon(1.0_dp)))
    allocate(w(n),oldbeta(size(beta)))
    do it=1,mit
      oldbeta=beta
      do i=1,n
        u=res(i)/max(s,tiny(1.0_dp))
        select case(trim(psifun))
        case("bisquare")
          if(abs(u)>tiny(1.0_dp)) then; w(i)=psi_bisquare(u)/u; else; w(i)=1.0_dp; end if
        case("hampel")
          if(abs(u)>tiny(1.0_dp)) then; w(i)=psi_hampel(u)/u; else; w(i)=1.0_dp; end if
        case default
          if(abs(u)>tiny(1.0_dp)) then; w(i)=psi_huber(u)/u; else; w(i)=1.0_dp; end if
        end select
        w(i)=max(0.0_dp,w(i))*cw(i)
      end do
      call weighted_least_squares(x,y,w,beta,res,cov,rank,status)
      if(trim(scalem)=="Huber" .or. trim(scalem)=="proposal 2") then
        s=sqrt(sum(min(res**2,(1.345_dp*s)**2))/max(0.7101645483_dp*real(n-rank,dp),1.0_dp))
      else
        s=mad_mass(res)
      end if
      delta=sqrt(sum((beta-oldbeta)**2)/max(sum(oldbeta**2),1.0e-20_dp))
      if(delta<tol) exit
    end do
    result%coefficients=beta;result%residuals=res;result%fitted=y-res;result%weights=w
    result%sigma=s;result%covariance=s*s*cov;result%rank=rank;result%df_residual=n-rank;result%iterations=it
    result%status=merge(mass_success,mass_no_convergence,it<=mit);result%method="robust linear model "//trim(meth)
    allocate(result%leverages(n));result%leverages=0.0_dp
  end subroutine rlm_fit

  function standardized_residuals(result) result(value)
    type(regression_result),intent(in)::result
    real(dp),allocatable::value(:)
    allocate(value(size(result%residuals)))
    value=result%residuals/(max(result%sigma,tiny(1.0_dp))*sqrt(max(1.0_dp-result%leverages,tiny(1.0_dp))))
  end function standardized_residuals

  function studentized_residuals(result) result(value)
    type(regression_result),intent(in)::result
    real(dp),allocatable::value(:),sr(:)
    real(dp)::df
    sr=standardized_residuals(result);allocate(value(size(sr)));df=real(result%df_residual,dp)
    value=sr/sqrt(max((df-sr*sr)/max(df-1.0_dp,1.0_dp),tiny(1.0_dp)))
  end function studentized_residuals

  subroutine boxcox_profile(x,y,lambdas,log_likelihood,status,eps)
    real(dp),intent(in)::x(:,:),y(:),lambdas(:)
    real(dp),allocatable,intent(out)::log_likelihood(:)
    integer,intent(out)::status
    real(dp),intent(in),optional::eps
    real(dp),allocatable::ys(:),yt(:),beta(:),res(:)
    real(dp)::e,la,gmean
    integer::i,rank,st,n
    n=size(y);e=0.02_dp;if(present(eps))e=eps
    if(size(x,1)/=n .or. any(y<=0.0_dp)) then;allocate(log_likelihood(0));status=mass_invalid_argument;return;end if
    gmean=exp(sum(log(y))/real(n,dp));ys=y/gmean;allocate(log_likelihood(size(lambdas)),yt(n))
    do i=1,size(lambdas)
      la=lambdas(i)
      if(abs(la)>e) then;yt=(ys**la-1.0_dp)/la
      else;yt=log(ys)*(1.0_dp+0.5_dp*la*log(ys)*(1.0_dp+la*log(ys)/3.0_dp*(1.0_dp+la*log(ys)/4.0_dp)));end if
      call least_squares(x,yt,beta,res,rank,st)
      log_likelihood(i)=-0.5_dp*real(n,dp)*log(max(sum(res**2),tiny(1.0_dp)))
    end do
    status=mass_success
  end subroutine boxcox_profile

  subroutine logtrans_profile(x,y,alpha,log_likelihood,status)
    real(dp),intent(in)::x(:,:),y(:),alpha(:)
    real(dp),allocatable,intent(out)::log_likelihood(:)
    integer,intent(out)::status
    real(dp),allocatable::yt(:),beta(:),res(:)
    integer::i,rank,st,n
    n=size(y)
    if(size(x,1)/=n .or. any(y+minval(alpha)<=0.0_dp)) then;allocate(log_likelihood(0));status=mass_invalid_argument;return;end if
    allocate(log_likelihood(size(alpha)),yt(n))
    do i=1,size(alpha)
      yt=log(y+alpha(i));call least_squares(x,yt,beta,res,rank,st)
      log_likelihood(i)=-0.5_dp*real(n,dp)*log(max(sum(res**2),tiny(1.0_dp)))-sum(yt)
    end do
    status=mass_success
  end subroutine logtrans_profile

  subroutine dose_p(coefficients,covariance,p,dose,se,status,link)
    real(dp),intent(in)::coefficients(:),covariance(:,:),p(:)
    real(dp),allocatable,intent(out)::dose(:),se(:)
    integer,intent(out)::status
    character(len=*),intent(in),optional::link
    real(dp)::eta,pd(2)
    character(len=16)::lnk
    integer::i
    lnk="logit";if(present(link))lnk=link
    if(size(coefficients)<2 .or. size(covariance,1)<2 .or. size(covariance,2)<2) then
      allocate(dose(0),se(0));status=mass_dimension_error;return
    end if
    allocate(dose(size(p)),se(size(p)))
    do i=1,size(p)
      select case(trim(lnk))
      case("probit");eta=inverse_normal_local(p(i))
      case("cloglog");eta=log(-log(max(1.0_dp-p(i),tiny(1.0_dp))))
      case default;eta=log(p(i)/max(1.0_dp-p(i),tiny(1.0_dp)))
      end select
      dose(i)=(eta-coefficients(1))/coefficients(2)
      pd=-[1.0_dp,dose(i)]/coefficients(2)
      se(i)=sqrt(max(0.0_dp,dot_product(pd,matmul(covariance(1:2,1:2),pd))))
    end do
    status=mass_success
  contains
    function inverse_normal_local(q) result(z)
      use mass_math,only:normal_quantile
      real(dp),intent(in)::q;real(dp)::z;z=normal_quantile(q)
    end function inverse_normal_local
  end subroutine dose_p

  subroutine huber_location(y,result,k,tolerance)
    real(dp),intent(in)::y(:)
    type(huber_result),intent(out)::result
    real(dp),intent(in),optional::k,tolerance
    real(dp)::kk,tol,mu,mu1,s
    integer::it
    kk=1.5_dp;if(present(k))kk=k;tol=1.0e-6_dp;if(present(tolerance))tol=tolerance
    if(size(y)<1)then;result%status=mass_invalid_argument;return;end if
    mu=median_local(y);s=mad_mass(y)
    if(s<=tiny(1.0_dp))then;result%status=mass_invalid_argument;return;end if
    do it=1,100
      mu1=sum(min(max(y,mu-kk*s),mu+kk*s))/real(size(y),dp)
      if(abs(mu-mu1)<tol*s)exit
      mu=mu1
    end do
    result%location=mu1;result%scale=s;result%iterations=it;result%status=mass_success
  contains
    function median_local(a) result(v)
      use mass_math,only:median_mass
      real(dp),intent(in)::a(:);real(dp)::v;v=median_mass(a)
    end function median_local
  end subroutine huber_location

  subroutine huber_location_scale(y,result,k,tolerance,fix_location,fix_scale)
    real(dp),intent(in)::y(:)
    type(huber_result),intent(out)::result
    real(dp),intent(in),optional::k,tolerance,fix_location,fix_scale
    real(dp)::kk,tol,mu,mu1,s,s1,th,beta
    logical::estmu,ests
    integer::it,n,n1
    kk=1.5_dp;if(present(k))kk=k;tol=1.0e-6_dp;if(present(tolerance))tol=tolerance;n=size(y)
    if(n<2)then;result%status=mass_invalid_argument;return;end if
    estmu=.not.present(fix_location);ests=.not.present(fix_scale)
    if(estmu)then;mu=median_internal(y);n1=n-1;else;mu=fix_location;n1=n;end if
    if(ests)then;s=mad_mass(y);else;s=fix_scale;end if
    if(s<=tiny(1.0_dp))then;result%location=mu;result%scale=0.0_dp;result%status=mass_success;return;end if
    th=2.0_dp*normal_cdf(kk)-1.0_dp;beta=th+kk*kk*(1.0_dp-th)-2.0_dp*kk*normal_pdf(kk)
    do it=1,30
      mu1=mu;if(estmu)mu1=sum(min(max(y,mu-kk*s),mu+kk*s))/real(n,dp)
      s1=s;if(ests)s1=sqrt(sum((min(max(y,mu-kk*s),mu+kk*s)-mu1)**2)/real(n1,dp)/beta)
      if(abs(mu-mu1)<tol*s .and. abs(s-s1)<tol*s)exit
      mu=mu1;s=s1
    end do
    result%location=mu1;result%scale=s1;result%iterations=it;result%status=mass_success
  contains
    function median_internal(a) result(v)
      use mass_math,only:median_mass
      real(dp),intent(in)::a(:);real(dp)::v;v=median_mass(a)
    end function median_internal
  end subroutine huber_location_scale

end module mass_regression
