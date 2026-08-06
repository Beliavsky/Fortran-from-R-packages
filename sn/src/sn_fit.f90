! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use sn_kinds, only : dp, tiny_dp
  use sn_status, only : sn_ok, sn_invalid_argument, sn_dimension_mismatch, &
                        sn_no_convergence
  use sn_linalg, only : inverse_general, sample_mean_covariance, covariance_to_correlation
  use sn_optim, only : nelder_mead, numerical_hessian
  use sn_math, only : normal_logpdf, normal_cdf, normal_quantile
  use sn_univariate, only : dsn, psn, dst, pst, dsc, psc, fournum, fournum_result
  use sn_multivariate, only : sn_mv_params, st_mv_params
  use sn_misc, only : q_penalty, galton_moors_to_alpha_nu
  implicit none
  private

  type, public :: selm_result
    character(len=8) :: family = 'SN'
    real(dp), allocatable :: beta(:)
    real(dp) :: omega = 1.0_dp
    real(dp) :: alpha = 0.0_dp
    real(dp) :: tau = 0.0_dp
    real(dp) :: nu = huge(1.0_dp)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: objective = huge(1.0_dp)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residuals(:)
    integer :: iterations = 0
    integer :: status = sn_ok
    logical :: converged = .false.
  end type selm_result

  type, public :: grouped_fit_result
    character(len=8) :: family = 'SN'
    real(dp) :: xi = 0.0_dp
    real(dp) :: omega = 1.0_dp
    real(dp) :: alpha = 0.0_dp
    real(dp) :: tau = 0.0_dp
    real(dp) :: nu = huge(1.0_dp)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp), allocatable :: covariance(:,:)
    integer :: iterations = 0
    integer :: status = sn_ok
    logical :: converged = .false.
  end type grouped_fit_result

  public :: selm_fit, predict_selm, fit_grouped
  public :: sn_mple, st_mple, msn_mle, msn_mple, mst_mple

contains

  subroutine selm_fit(x, y, family, fit, weights, penalty, max_iter, tol, start)
    real(dp), intent(in) :: x(:,:), y(:)
    character(len=*), intent(in), optional :: family
    type(selm_result), intent(out) :: fit
    real(dp), intent(in), optional :: weights(:)
    character(len=*), intent(in), optional :: penalty
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol, start(:)
    real(dp), allocatable :: w(:), theta0(:), theta(:), xtwx(:,:), xtwy(:), inv(:,:)
    real(dp), allocatable :: hess(:,:), cov(:,:), beta0(:), r(:)
    real(dp) :: fbest, scale0, skew0, sumw, opt_tol
    integer :: n, p, k, ierr, iters, maxit, npar, i
    character(len=8) :: fam
    character(len=8) :: pen

    n = size(x,1)
    p = size(x,2)
    fam = 'SN'
    if (present(family)) fam = uppercase(adjustl(family))
    pen = 'NONE'
    if (present(penalty)) pen = uppercase(adjustl(penalty))
    fit%family = fam
    if (n < 2 .or. p < 1 .or. size(y) /= n) then
      fit%status = sn_dimension_mismatch
      return
    end if
    allocate(w(n))
    w = 1.0_dp
    if (present(weights)) then
      if (size(weights) /= n .or. any(weights < 0.0_dp)) then
        fit%status = sn_invalid_argument
        return
      end if
      w = weights
    end if
    sumw = sum(w)
    if (sumw <= 0.0_dp) then
      fit%status = sn_invalid_argument
      return
    end if
    allocate(xtwx(p,p),xtwy(p),beta0(p),r(n))
    xtwx = 0.0_dp
    xtwy = 0.0_dp
    do i=1,n
      xtwx = xtwx+w(i)*outer(x(i,:),x(i,:))
      xtwy = xtwy+w(i)*x(i,:)*y(i)
    end do
    call inverse_general(xtwx,inv,ierr)
    if (ierr == sn_ok) then
      beta0 = matmul(inv,xtwy)
    else
      beta0 = 0.0_dp
      if (p >= 1) beta0(1) = sum(w*y)/sumw
    end if
    r = y-matmul(x,beta0)
    scale0 = sqrt(max(tiny_dp,sum(w*r*r)/sumw))
    skew0 = weighted_skewness(r,w)

    select case (trim(fam))
    case ('NORMAL','N')
      npar = p+1
    case ('SN')
      npar = p+2
    case ('ESN')
      npar = p+3
    case ('ST')
      npar = p+3
    case ('SC')
      npar = p+2
    case default
      fit%status = sn_invalid_argument
      return
    end select
    allocate(theta0(npar))
    theta0(1:p) = beta0
    theta0(p+1) = log(scale0)
    if (npar > p+1) theta0(p+2) = sign(min(3.0_dp,4.0_dp*abs(skew0)),skew0)
    if (trim(fam) == 'ESN') theta0(p+3) = 0.0_dp
    if (trim(fam) == 'ST') theta0(p+3) = log(8.0_dp)
    if (present(start)) then
      if (size(start) == npar) theta0 = start
    end if
    maxit = 5000
    if (present(max_iter)) maxit = max_iter
    opt_tol = 1.0e-8_dp
    if (present(tol)) opt_tol = tol
    call nelder_mead(objective,theta0,theta,fbest,ierr,iters,max_iter=maxit, &
                     tol=opt_tol,step=0.15_dp)
    call unpack(theta,fit)
    fit%objective = fbest
    fit%log_likelihood = loglik(theta)
    fit%iterations = iters
    fit%status = ierr
    fit%converged = ierr == sn_ok
    allocate(fit%fitted(n),fit%residuals(n))
    fit%fitted = matmul(x,fit%beta)
    fit%residuals = y-fit%fitted
    call numerical_hessian(objective,theta,hess,cov,k,step=2.0e-4_dp)
    if (k == sn_ok) then
      fit%covariance = cov
    else
      allocate(fit%covariance(0,0))
    end if

  contains
    real(dp) function objective(t) result(value)
      real(dp), intent(in) :: t(:)
      value = -loglik(t)
      if (trim(pen) == 'Q' .and. size(t) >= p+2) then
        if (trim(fam) == 'ST') then
          value = value+q_penalty(t(p+2),1.0_dp+exp(t(p+3)))
        else
          value = value+q_penalty(t(p+2))
        end if
      end if
      if (.not. ieee_is_finite(value)) value = 0.25_dp*huge(1.0_dp)
    end function objective

    real(dp) function loglik(t) result(value)
      real(dp), intent(in) :: t(:)
      real(dp) :: z, om, a, ta, df, lp
      integer :: j
      om = exp(t(p+1))
      a = 0.0_dp
      ta = 0.0_dp
      df = huge(1.0_dp)
      if (size(t)>=p+2) a = t(p+2)
      if (trim(fam)=='ESN') ta = t(p+3)
      if (trim(fam)=='ST') df = 1.0_dp+exp(t(p+3))
      value = 0.0_dp
      do j=1,n
        z = y(j)-dot_product(x(j,:),t(1:p))
        lp = -huge(1.0_dp)
        select case (trim(fam))
        case ('NORMAL','N')
          lp = normal_logpdf(z/om)-log(om)
        case ('SN')
          lp = dsn(z,omega=om,alpha=a,log_pdf=.true.)
        case ('ESN')
          lp = dsn(z,omega=om,alpha=a,tau=ta,log_pdf=.true.)
        case ('ST')
          lp = dst(z,omega=om,alpha=a,nu=df,log_pdf=.true.)
        case ('SC')
          lp = dsc(z,omega=om,alpha=a,log_pdf=.true.)
        case default
          lp = -huge(1.0_dp)
        end select
        if (.not. ieee_is_finite(lp)) then
          value = -huge(1.0_dp)
          return
        end if
        value = value+w(j)*lp
      end do
    end function loglik

    subroutine unpack(t,out)
      real(dp), intent(in) :: t(:)
      type(selm_result), intent(inout) :: out
      allocate(out%beta(p))
      out%beta = t(1:p)
      out%omega = exp(t(p+1))
      out%alpha = 0.0_dp
      out%tau = 0.0_dp
      out%nu = huge(1.0_dp)
      if (size(t)>=p+2) out%alpha = t(p+2)
      if (trim(fam)=='ESN') out%tau = t(p+3)
      if (trim(fam)=='ST') out%nu = 1.0_dp+exp(t(p+3))
      if (trim(fam)=='SC') out%nu = 1.0_dp
    end subroutine unpack
  end subroutine selm_fit

  subroutine predict_selm(fit, x, location, mean_response, info)
    type(selm_result), intent(in) :: fit
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: location(:)
    real(dp), allocatable, intent(out), optional :: mean_response(:)
    integer, intent(out), optional :: info
    real(dp) :: delta, bnu
    if (.not. allocated(fit%beta) .or. size(x,2)/=size(fit%beta)) then
      allocate(location(0))
      if (present(mean_response)) allocate(mean_response(0))
      if (present(info)) info=sn_dimension_mismatch
      return
    end if
    allocate(location(size(x,1)))
    location=matmul(x,fit%beta)
    if (present(mean_response)) then
      allocate(mean_response(size(location)))
      mean_response=location
      delta=fit%alpha/sqrt(1.0_dp+fit%alpha*fit%alpha)
      select case(trim(fit%family))
      case('SN')
        mean_response=mean_response+fit%omega*delta*sqrt(2.0_dp/acos(-1.0_dp))
      case('ST')
        if (fit%nu>1.0_dp) then
          bnu=sqrt(fit%nu/acos(-1.0_dp))*exp(log_gamma(0.5_dp*(fit%nu-1.0_dp))- &
                log_gamma(0.5_dp*fit%nu))
          mean_response=mean_response+fit%omega*delta*bnu
        end if
      end select
    end if
    if (present(info)) info=sn_ok
  end subroutine predict_selm

  subroutine sn_mple(y, fit, max_iter, tol, penalty)
    real(dp), intent(in) :: y(:)
    type(selm_result), intent(out) :: fit
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    character(len=*), intent(in), optional :: penalty
    real(dp), allocatable :: x(:,:)
    integer :: mi
    real(dp) :: tt
    character(len=8) :: pp
    allocate(x(size(y),1)); x(:,1)=1.0_dp
    mi=5000; if(present(max_iter)) mi=max_iter
    tt=1.0e-8_dp; if(present(tol)) tt=tol
    pp='Q'; if(present(penalty)) pp=penalty
    call selm_fit(x,y,'SN',fit,penalty=pp,max_iter=mi,tol=tt)
  end subroutine sn_mple

  subroutine st_mple(y, fit, max_iter, tol, penalty)
    real(dp), intent(in) :: y(:)
    type(selm_result), intent(out) :: fit
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    character(len=*), intent(in), optional :: penalty
    real(dp), allocatable :: x(:,:)
    integer :: mi
    real(dp) :: tt
    character(len=8) :: pp
    allocate(x(size(y),1)); x(:,1)=1.0_dp
    mi=5000; if(present(max_iter)) mi=max_iter
    tt=1.0e-8_dp; if(present(tol)) tt=tol
    pp='Q'; if(present(penalty)) pp=penalty
    call selm_fit(x,y,'ST',fit,penalty=pp,max_iter=mi,tol=tt)
  end subroutine st_mple

  subroutine msn_mle(data, params, info, fits)
    real(dp), intent(in) :: data(:,:)
    type(sn_mv_params), intent(out) :: params
    integer, intent(out) :: info
    type(selm_result), allocatable, intent(out), optional :: fits(:)
    call composite_mv_fit(data,.false.,params=params,info=info,snfits=fits)
  end subroutine msn_mle

  subroutine msn_mple(data, params, info, fits)
    real(dp), intent(in) :: data(:,:)
    type(sn_mv_params), intent(out) :: params
    integer, intent(out) :: info
    type(selm_result), allocatable, intent(out), optional :: fits(:)
    call composite_mv_fit(data,.false.,params=params,info=info,snfits=fits,penalized=.true.)
  end subroutine msn_mple

  subroutine mst_mple(data, params, info, fits)
    real(dp), intent(in) :: data(:,:)
    type(st_mv_params), intent(out) :: params
    integer, intent(out) :: info
    type(selm_result), allocatable, intent(out), optional :: fits(:)
    call composite_mv_fit(data,.true.,stparams=params,info=info,snfits=fits,penalized=.true.)
  end subroutine mst_mple

  subroutine composite_mv_fit(data,is_t,params,stparams,info,snfits,penalized)
    real(dp), intent(in) :: data(:,:)
    logical, intent(in) :: is_t
    type(sn_mv_params), intent(out), optional :: params
    type(st_mv_params), intent(out), optional :: stparams
    integer, intent(out) :: info
    type(selm_result), allocatable, intent(out), optional :: snfits(:)
    logical, intent(in), optional :: penalized
    type(selm_result), allocatable :: f(:)
    real(dp), allocatable :: mean(:),cov(:,:),cor(:,:),sd(:),omega(:,:),xi(:),alpha(:),nus(:)
    real(dp) :: delta, shift
    integer :: d,j,ierr
    character(len=8) :: pen
    d=size(data,2)
    if (size(data,1)<3 .or. d<1) then
      info=sn_invalid_argument
      return
    end if
    pen='NONE'
    if (present(penalized)) then
      if (penalized) pen='Q'
    end if
    allocate(f(d),xi(d),alpha(d),nus(d))
    do j=1,d
      if (is_t) then
        call fit_intercept(data(:,j),'ST',f(j),pen)
      else
        call fit_intercept(data(:,j),'SN',f(j),pen)
      end if
      if (.not.f(j)%converged) then
        info=sn_no_convergence
        return
      end if
      delta=f(j)%alpha/sqrt(1.0_dp+f(j)%alpha*f(j)%alpha)
      if (is_t .and. f(j)%nu>1.0_dp) then
        shift=f(j)%omega*delta*sqrt(f(j)%nu/acos(-1.0_dp))* &
          exp(log_gamma(0.5_dp*(f(j)%nu-1.0_dp))-log_gamma(0.5_dp*f(j)%nu))
      else
        shift=f(j)%omega*delta*sqrt(2.0_dp/acos(-1.0_dp))
      end if
      xi(j)=f(j)%beta(1)
      alpha(j)=f(j)%alpha
      nus(j)=f(j)%nu
    end do
    call sample_mean_covariance(data,mean,cov,ierr)
    call covariance_to_correlation(cov,cor,sd,ierr)
    if (ierr/=sn_ok) then
      info=ierr
      return
    end if
    allocate(omega(d,d))
    do j=1,d
      do ierr=1,d
        omega(ierr,j)=cor(ierr,j)*f(ierr)%omega*f(j)%omega
      end do
    end do
    if (is_t) then
      if (.not.present(stparams)) then
        info=sn_invalid_argument
        return
      end if
      stparams%xi=xi; stparams%omega=omega; stparams%alpha=alpha
      stparams%nu=sum(nus)/real(d,dp)
    else
      if (.not.present(params)) then
        info=sn_invalid_argument
        return
      end if
      params%xi=xi; params%omega=omega; params%alpha=alpha; params%tau=0.0_dp
    end if
    if (present(snfits)) snfits=f
    info=sn_ok
  contains
    subroutine fit_intercept(v,family0,out,pen0)
      real(dp), intent(in) :: v(:)
      character(len=*), intent(in) :: family0,pen0
      type(selm_result), intent(out) :: out
      real(dp), allocatable :: xx(:,:)
      allocate(xx(size(v),1)); xx(:,1)=1.0_dp
      call selm_fit(xx,v,family0,out,penalty=pen0,max_iter=3500,tol=2.0e-7_dp)
    end subroutine fit_intercept
  end subroutine composite_mv_fit

  subroutine fit_grouped(lower, upper, counts, family, fit, max_iter, tol)
    real(dp), intent(in) :: lower(:),upper(:),counts(:)
    character(len=*), intent(in), optional :: family
    type(grouped_fit_result), intent(out) :: fit
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: theta0(:),theta(:),h(:,:),cov(:,:)
    real(dp) :: fbest,midmean,midsd,total,opt_tol
    integer :: m,npar,ierr,iters,i,maxit
    character(len=8) :: fam
    fam='SN'; if(present(family)) fam=uppercase(adjustl(family))
    m=size(lower)
    fit%family=fam
    if(size(upper)/=m .or. size(counts)/=m .or. m<1 .or. any(upper<=lower) .or. any(counts<0.0_dp)) then
      fit%status=sn_invalid_argument; return
    end if
    total=sum(counts)
    if(total<=0.0_dp) then
      fit%status=sn_invalid_argument; return
    end if
    midmean=sum(counts*0.5_dp*(lower+upper))/total
    midsd=sqrt(max(tiny_dp,sum(counts*(0.5_dp*(lower+upper)-midmean)**2)/total))
    select case(trim(fam))
    case('NORMAL','N'); npar=2
    case('SN','SC'); npar=3
    case('ESN','ST'); npar=4
    case default; fit%status=sn_invalid_argument; return
    end select
    allocate(theta0(npar)); theta0=0.0_dp
    theta0(1)=midmean; theta0(2)=log(midsd)
    if(trim(fam)=='ST') theta0(4)=log(8.0_dp)
    maxit=4000; if(present(max_iter)) maxit=max_iter
    opt_tol=1.0e-8_dp; if(present(tol)) opt_tol=tol
    call nelder_mead(obj,theta0,theta,fbest,ierr,iters,max_iter=maxit,tol=opt_tol,step=0.15_dp)
    fit%xi=theta(1); fit%omega=exp(theta(2)); fit%alpha=0.0_dp; fit%tau=0.0_dp; fit%nu=huge(1.0_dp)
    if(npar>=3) fit%alpha=theta(3)
    if(trim(fam)=='ESN') fit%tau=theta(4)
    if(trim(fam)=='ST') fit%nu=1.0_dp+exp(theta(4))
    if(trim(fam)=='SC') fit%nu=1.0_dp
    fit%log_likelihood=-fbest; fit%iterations=iters; fit%status=ierr; fit%converged=ierr==sn_ok
    call numerical_hessian(obj,theta,h,cov,i,step=2.0e-4_dp)
    if(i==sn_ok) then
      fit%covariance=cov
    else
      allocate(fit%covariance(0,0))
    end if
  contains
    real(dp) function obj(t) result(v)
      real(dp),intent(in)::t(:)
      real(dp)::pl,pu,prob,om,a,ta,df
      integer::j
      om=exp(t(2)); a=0.0_dp; ta=0.0_dp; df=huge(1.0_dp)
      if(size(t)>=3) a=t(3)
      if(trim(fam)=='ESN') ta=t(4)
      if(trim(fam)=='ST') df=1.0_dp+exp(t(4))
      v=0.0_dp
      do j=1,m
        pl=0.0_dp; pu=0.0_dp
        select case(trim(fam))
        case('NORMAL','N')
          pu=normal_cdf((upper(j)-t(1))/om); pl=normal_cdf((lower(j)-t(1))/om)
        case('SN')
          pu=psn(upper(j),t(1),om,a); pl=psn(lower(j),t(1),om,a)
        case('ESN')
          pu=psn(upper(j),t(1),om,a,ta); pl=psn(lower(j),t(1),om,a,ta)
        case('ST')
          pu=pst(upper(j),t(1),om,a,df); pl=pst(lower(j),t(1),om,a,df)
        case('SC')
          pu=psc(upper(j),t(1),om,a); pl=psc(lower(j),t(1),om,a)
        case default
          pu=0.0_dp; pl=0.0_dp
        end select
        prob=max(tiny_dp,pu-pl)
        v=v-counts(j)*log(prob)
      end do
      if(.not.ieee_is_finite(v)) v=0.25_dp*huge(1.0_dp)
    end function obj
  end subroutine fit_grouped

  pure function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:)
    real(dp)::c(size(a),size(b))
    integer::i,j
    do j=1,size(b); do i=1,size(a); c(i,j)=a(i)*b(j); end do; end do
  end function outer

  pure real(dp) function weighted_skewness(x,w) result(s)
    real(dp),intent(in)::x(:),w(:)
    real(dp)::mu,v,total
    total=sum(w)
    mu=sum(w*x)/max(total,tiny_dp)
    v=sum(w*(x-mu)**2)/max(total,tiny_dp)
    if(v<=tiny_dp) then
      s=0.0_dp
    else
      s=sum(w*(x-mu)**3)/max(total,tiny_dp)/v**1.5_dp
    end if
  end function weighted_skewness

  pure function uppercase(text) result(out)
    character(len=*),intent(in)::text
    character(len=len(text))::out
    integer::i,k
    out=text
    do i=1,len(text)
      k=iachar(out(i:i))
      if(k>=iachar('a').and.k<=iachar('z')) out(i:i)=achar(k-32)
    end do
  end function uppercase

end module sn_fit
