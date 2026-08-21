! SPDX-License-Identifier: GPL-2.0-or-later
module tsa_arma
  use tsa_kinds, only : dp
  use tsa_types, only : ar_fit_result, spectrum_result
  use tsa_utils, only : mean_value, variance_n, build_lag_matrix, linear_fit, sort_with_index
  use tsa_statistics, only : autocorrelation
  use tseries_linalg, only : least_squares, solve_linear
  use leaps, only : regsubsets_result, regsubsets_fit
  implicit none
  private
  public :: ar_ols_fit, arma_spectrum, eacf, armasubsets_fit, boxcox_ar
  public :: prewhiten_filter

contains

  function ar_ols_fit(x,order_max,select_aic,intercept) result(fit)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: order_max
    logical, intent(in), optional :: select_aic,intercept
    type(ar_fit_result) :: fit
    type(ar_fit_result) :: trial
    logical :: sel,inc
    integer :: p,maxp
    real(dp) :: best
    sel=.true.
    if(present(select_aic))sel=select_aic
    inc=.true.
    if(present(intercept))inc=intercept
    maxp=min(order_max,max(0,size(x)-3))
    best=huge(1.0_dp)
    if(.not.sel) then
      call fit_order(maxp,fit)
      return
    end if
    do p=0,maxp
      call fit_order(p,trial)
      if(trial%status==0 .and. trial%aic<best) then
        best=trial%aic
        fit=trial
      end if
    end do
  contains
    subroutine fit_order(p,r)
      integer,intent(in)::p
      type(ar_fit_result),intent(out)::r
      real(dp),allocatable::y(:),xx(:,:),beta(:),resid(:)
      integer::status,nobs,k
      r%order=p
      r%intercept=inc
      if(p==0) then
        allocate(r%coefficients(merge(1,0,inc)),r%residuals(size(x)),r%fitted(size(x)))
        if(inc)then
        r%coefficients(1)=mean_value(x)
        r%fitted=r%coefficients(1)
        else
        r%fitted=0.0_dp
        end if
        r%residuals=x-r%fitted
        r%variance=sum(r%residuals**2)/real(size(x),dp)
        k=merge(1,0,inc)
        r%aic=real(size(x),dp)*log(max(r%variance,tiny(1.0_dp)))+2.0_dp*real(k,dp)
        r%status=0
        return
      end if
      call build_lag_matrix(x,p,p+1,y,xx,inc)
      allocate(beta(size(xx,2)),resid(size(y)))
      call least_squares(xx,y,beta,residuals=resid,status=status)
      r%status=status
      if(status/=0)return
      allocate(r%coefficients(size(beta)),r%residuals(size(x)),r%fitted(size(x)))
      r%coefficients=beta
      r%residuals=0.0_dp
      r%fitted=x
      r%residuals(p+1:)=resid
      r%fitted(p+1:)=y-resid
      nobs=size(y)
      r%variance=sum(resid**2)/real(nobs,dp)
      k=size(beta)
      r%aic=real(nobs,dp)*log(max(r%variance,tiny(1.0_dp)))+2.0_dp*real(k,dp)
    end subroutine fit_order
  end function ar_ols_fit

  function arma_spectrum(freq,ar,ma,sigma2,sar,sma,period) result(res)
    real(dp), intent(in) :: freq(:)
    real(dp), intent(in), optional :: ar(:),ma(:),sigma2,sar(:),sma(:)
    integer, intent(in), optional :: period
    type(spectrum_result) :: res
    integer :: s
    real(dp)::v
    allocate(res%frequency(size(freq)), res%spectrum(size(freq)))
    res%frequency = freq
    res%spectrum = 1.0_dp
    s = 12
    if (present(period)) s = period
    if(present(ar))res%spectrum=res%spectrum/poly_power(freq,ar,1,.true.)
    if(present(ma))res%spectrum=res%spectrum*poly_power(freq,ma,1,.false.)
    if(present(sar))res%spectrum=res%spectrum/poly_power(freq,sar,s,.true.)
    if(present(sma))res%spectrum=res%spectrum*poly_power(freq,sma,s,.false.)
    v=1.0_dp
    if(present(sigma2))v=sigma2
    res%spectrum=res%spectrum*v
  contains
    function poly_power(f,coef,per,is_ar) result(pow)
      real(dp),intent(in)::f(:),coef(:)
      integer,intent(in)::per
      logical,intent(in)::is_ar
      real(dp)::pow(size(f)),re,im,ang,sgn
      integer::ii,jj
      sgn=merge(-1.0_dp,1.0_dp,is_ar)
      do ii=1,size(f)
        re=1.0_dp
        im=0.0_dp
        do jj=1,size(coef)
          ang=2.0_dp*acos(-1.0_dp)*f(ii)*real(jj*per,dp)
          re=re+sgn*coef(jj)*cos(ang)
          im=im-sgn*coef(jj)*sin(ang)
        end do
        pow(ii)=re*re+im*im
      end do
    end function poly_power
  end function arma_spectrum

  subroutine prewhiten_filter(x,y,ar,x_filtered,y_filtered)
    real(dp),intent(in)::x(:),y(:),ar(:)
    real(dp),allocatable,intent(out)::x_filtered(:),y_filtered(:)
    integer::n,i,j,p
    n=min(size(x),size(y))
    p=size(ar)
    allocate(x_filtered(n),y_filtered(n))
    x_filtered=0.0_dp
    y_filtered=0.0_dp
    do i=1,n
      x_filtered(i)=x(i)
      y_filtered(i)=y(i)
      do j=1,min(p,i-1)
        x_filtered(i)=x_filtered(i)-ar(j)*x(i-j)
        y_filtered(i)=y_filtered(i)-ar(j)*y(i-j)
      end do
    end do
  end subroutine prewhiten_filter

  subroutine armasubsets_fit(y,nar,nma,result,nvmax,nbest,method,status)
    real(dp),intent(in)::y(:)
    integer,intent(in)::nar,nma
    type(regsubsets_result),intent(out)::result
    integer,intent(in),optional::nvmax,nbest
    character(len=*),intent(in),optional::method
    integer,intent(out),optional::status
    type(ar_fit_result)::af
    real(dp),allocatable::xx(:,:),yy(:)
    integer::i,j,start,n,ier,nv,nb
    character(len=16)::meth
    af=ar_ols_fit(y,min(max(1,nar),max(1,int(sqrt(real(size(y),dp))))),select_aic=.true.)
    start=max(nar,nma)+1
    n=size(y)-start+1
    if(n<=2 .or. nar+nma<1)then
    ier=1
    if(present(status))status=ier
    return
    end if
    allocate(xx(n,nar+nma),yy(n))
    yy=y(start:)
    do i=1,n
    do j=1,nar
    xx(i,j)=y(start+i-1-j)
    end do
    do j=1,nma
    xx(i,nar+j)=af%residuals(start+i-1-j)
    end do
    end do
    nv=min(8,nar+nma)
    if(present(nvmax))nv=nvmax
    nb=1
    if(present(nbest))nb=nbest
    meth='exhaustive'
    if(present(method))meth=method
    call regsubsets_fit(xx,yy,result,nvmax=nv,nbest=nb,method=trim(meth),intercept=.true.,ier=ier)
    if(present(status))status=ier
  end subroutine armasubsets_fit

  subroutine boxcox_ar(y,lambda,order,loglik,mle,ci,status)
    real(dp),intent(in)::y(:),lambda(:)
    integer,intent(in),optional::order
    real(dp),allocatable,intent(out)::loglik(:)
    real(dp),intent(out)::mle,ci(2)
    integer,intent(out)::status
    real(dp),allocatable::ys(:)
    type(ar_fit_result)::af
    integer::i,ord,imax
    real(dp)::scale,sumlog,limit
    status=0
    if(any(y<=0.0_dp))then
    status=1
    allocate(loglik(0))
    mle=0.0_dp
    ci=0.0_dp
    return
    end if
    scale=maxval(abs(y))+1.0_dp
    allocate(loglik(size(lambda)))
    sumlog=sum(log(y/scale))
    ord=-1
    if(present(order))ord=order
    do i=1,size(lambda)
      allocate(ys(size(y)))
      if(abs(lambda(i))>epsilon(1.0_dp))then
      ys=((y/scale)**lambda(i)-1.0_dp)/lambda(i)
      else
      ys=log(y/scale)
      end if
      if(ord>=0)then
      af=ar_ols_fit(ys,ord,select_aic=.false.)
      else
      af=ar_ols_fit(ys,min(20,size(y)/5),select_aic=.true.)
      end if
      loglik(i)=-0.5_dp*real(size(y),dp)*log(max(af%variance,tiny(1.0_dp)))+(lambda(i)-1.0_dp)*sumlog
      deallocate(ys)
    end do
    imax=maxloc(loglik,dim=1)
    mle=lambda(imax)
    limit=loglik(imax)-1.920729410347062_dp
    ci(1)=lambda(1)
    ci(2)=lambda(size(lambda))
    do i=1,size(lambda)
    if(loglik(i)>=limit)then
    ci(1)=lambda(i)
    exit
    end if
    end do
    do i=size(lambda),1,-1
    if(loglik(i)>=limit)then
    ci(2)=lambda(i)
    exit
    end if
    end do
  end subroutine boxcox_ar

  subroutine eacf(z,ar_max,ma_max,eacfm,symbol,status)
    real(dp),intent(in)::z(:)
    integer,intent(in)::ar_max,ma_max
    real(dp),allocatable,intent(out)::eacfm(:,:)
    character(len=1),allocatable,intent(out)::symbol(:,:)
    integer,intent(out)::status
    real(dp),allocatable::zc(:),zm(:,:),covp(:),cov1(:),m1(:,:),m2(:,:),temp(:),ac(:)
    type(ar_fit_result)::af
    integer::nar,nma,ncov,nrow,ncol,i,j,i1,count,n,ii
    real(dp)::den
    status=0
    nar=ar_max
    nma=ma_max+1
    ncov=nar+nma+2
    nrow=nar+nma+1
    ncol=nrow-1
    n=size(z)
    if(n<=nrow+2)then
    status=1
    allocate(eacfm(0,0),symbol(0,0))
    return
    end if
    allocate(zc(n))
    zc=z-mean_value(z)
    allocate(zm(n,nar))
    zm=0.0_dp
    do i=1,nar
    if(i<n)zm(i+1:,i)=zc(:n-i)
    end do
    call autocorrelation(zc,ncov,covp,drop_lag_zero=.false.)
    allocate(cov1(2*ncov+1))
    do i=-ncov,ncov
    cov1(i+ncov+1)=covp(abs(i)+1)
    end do
    allocate(m1(nrow,ncol))
    m1=0.0_dp
    do i=1,ncol
      af=ar_ols_fit(zc,i,select_aic=.false.,intercept=.false.)
      if(af%status/=0)then
      status=2
      return
      end if
      m1(1:i,i)=af%coefficients(1:i)
    end do
    allocate(eacfm(nar+1,ma_max+1))
    eacfm=0.0_dp
    do count=1,nma
      if(ncol<=1)exit
      allocate(m2(nrow,ncol-1))
      m2=0.0_dp
      do i=1,ncol-1
        i1=i+1
        allocate(temp(nrow))
        temp=0.0_dp
        temp(2:)=m1(:nrow-1,i)
        temp(1)=-1.0_dp
        den=m1(i,i)
        if(abs(den)<=tiny(1.0_dp))den=sign(tiny(1.0_dp),den+tiny(1.0_dp))
        m2(:,i)=m1(:,i1)-temp*m1(i1,i1)/den
        m2(i1,i)=0.0_dp
        deallocate(temp)
      end do
      eacfm(1,count)=cov1(ncov+1+count)
      do i=1,nar
        allocate(temp(n-i))
        temp=zc(i+1:)
        do j=1,i
        temp=temp-m2(j,i)*zm(i+1:,j)
        end do
        call autocorrelation(temp,count,ac,drop_lag_zero=.false.)
        eacfm(i+1,count)=ac(count+1)
        deallocate(temp,ac)
      end do
      call move_alloc(m2,m1)
      ncol=ncol-1
    end do
    allocate(symbol(nar+1,ma_max+1))
    do j=1,ma_max+1
    do i=1,nar+1
    ii=n-i+1-j
    symbol(i,j)=merge('x','o',abs(eacfm(i,j))>2.0_dp/sqrt(real(max(1,ii),dp)))
    end do
    end do
  end subroutine eacf
end module tsa_arma
