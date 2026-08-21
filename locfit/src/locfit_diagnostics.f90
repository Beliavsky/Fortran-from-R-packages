! Residuals and simple fit criteria derived from locfit R/locfit.r and src/fitted.c.
! GPL-2-or-later; see LICENSE_NOTICE and upstream/locfit-R.
module locfit_diagnostics
  use locfit_kinds, only : dp
  use locfit_constants
  use locfit_families, only : family_terms
  implicit none
  private
  public :: locfit_residual, studentize_residual
  public :: gcv_score, aic_score, cp_score
  public :: km_mean_residual_life, median_value
contains

  real(dp) function locfit_residual(y,w,theta,family,link,residual_type,status, &
      censored,robust_scale) result(value)
    real(dp),intent(in)::y,w,theta
    integer,intent(in)::family,link,residual_type
    integer,intent(out)::status
    logical,intent(in),optional::censored
    real(dp),intent(in),optional::robust_scale
    real(dp)::r(llen),raw,rs
    logical::ce
    ce=.false.;if(present(censored))ce=censored
    rs=1.0_dp;if(present(robust_scale))rs=robust_scale
    call family_terms(theta,y,family,link,r,status,censored=ce,prior_weight=w,robust_scale=rs)
    if(status/=lf_ok)then
      value=0.0_dp
      return
    end if
    select case(iand(family,63))
    case(tgaus,trobt,tcauc)
      raw=y-r(zmean)
    case default
      raw=y-w*r(zmean)
    end select
    select case(residual_type)
    case(rdev)
      value=sqrt(max(0.0_dp,-2.0_dp*r(zlik)))
      if(r(zdll)<=0.0_dp)value=-value
    case(rpear)
      if(r(zddll)<=0.0_dp)then
        value=merge(0.0_dp,huge(1.0_dp),r(zdll)==0.0_dp)
      else
        value=r(zdll)/sqrt(r(zddll))
      end if
    case(rraw)
      value=raw
    case(rldot)
      value=r(zdll)
    case(rdev2)
      value=-2.0_dp*r(zlik)
    case(rlddt)
      value=r(zddll)
    case(rfit)
      value=theta
    case(rmean)
      value=r(zmean)
    case default
      status=lf_err
      value=0.0_dp
    end select
  end function locfit_residual

  pure real(dp) function studentize_residual(res,influence,nlx,link_zdd,residual_type) result(value)
    real(dp),intent(in)::res,influence,nlx,link_zdd
    integer,intent(in)::residual_type
    real(dp)::inl,var,den
    inl=min(1.0_dp,influence*link_zdd)
    var=min(inl,nlx*nlx*link_zdd)
    den=1.0_dp-2.0_dp*inl+var
    if(den<0.0_dp)then
      value=0.0_dp
      return
    end if
    select case(residual_type)
    case(rdev,rpear,rraw,rldot)
      value=res/sqrt(max(den,tiny(1.0_dp)))
    case(rdev2)
      value=res/max(den,tiny(1.0_dp))
    case default
      value=res
    end select
  end function studentize_residual

  pure real(dp) function gcv_score(n,loglik,df1) result(value)
    integer,intent(in)::n
    real(dp),intent(in)::loglik,df1
    real(dp)::den
    den=real(n,dp)-df1
    if(abs(den)<=sqrt(epsilon(1.0_dp)))then
      value=huge(1.0_dp)
    else
      value=(-2.0_dp*real(n,dp)*loglik)/(den*den)
    end if
  end function gcv_score

  pure real(dp) function aic_score(loglik,df1,penalty) result(value)
    real(dp),intent(in)::loglik,df1
    real(dp),intent(in),optional::penalty
    real(dp)::pen
    pen=2.0_dp;if(present(penalty))pen=penalty
    value=-2.0_dp*loglik+pen*df1
  end function aic_score

  pure real(dp) function cp_score(n,loglik,df1,sigma2) result(value)
    integer,intent(in)::n
    real(dp),intent(in)::loglik,df1,sigma2
    if(sigma2<=0.0_dp)then
      value=huge(1.0_dp)
    else
      value=(-2.0_dp*loglik)/sigma2-real(n,dp)+2.0_dp*df1
    end if
  end function cp_score

  pure real(dp) function median_value(x) result(value)
    real(dp),intent(in)::x(:)
    real(dp),allocatable::a(:)
    real(dp)::key
    integer::i,j,n
    n=size(x)
    if(n==0)then;value=0.0_dp;return;end if
    allocate(a(n));a=x
    do i=2,n
      key=a(i);j=i-1
      do while(j>=1)
        if(a(j)<=key)exit
        a(j+1)=a(j);j=j-1
      end do
      a(j+1)=key
    end do
    if(mod(n,2)==1)then
      value=a((n+1)/2)
    else
      value=0.5_dp*(a(n/2)+a(n/2+1))
    end if
  end function median_value

  pure subroutine km_mean_residual_life(times,censored,mrl)
    ! Translation of R km.mrl(). censored=.true. denotes a censored observation.
    real(dp),intent(in)::times(:)
    logical,intent(in)::censored(:)
    real(dp),intent(out)::mrl(:)
    integer::n,i,j,keyi
    integer,allocatable::ord(:)
    real(dp),allocatable::ts(:),surv(:),ints(:),mr(:)
    logical,allocatable::cs(:)
    real(dp)::keys
    n=size(times);mrl=0.0_dp
    if(size(censored)/=n .or. size(mrl)/=n .or. n==0)return
    allocate(ord(n),ts(n),cs(n),surv(max(1,n-1)),ints(max(1,n-1)),mr(n))
    cs=.false.;surv=0.0_dp;ints=0.0_dp;mr=0.0_dp
    do i=1,n;ord(i)=i;end do
    do i=2,n
      keyi=ord(i);keys=times(keyi);j=i-1
      do while(j>=1)
        if(times(ord(j))<=keys)exit
        ord(j+1)=ord(j);j=j-1
      end do
      ord(j+1)=keyi
    end do
    do i=1,n
      ts(i)=times(ord(i));cs(i)=censored(ord(i))
    end do
    if(n==1)then
      mrl(1)=0.0_dp
      return
    end if
    surv(1)=1.0_dp-(merge(0.0_dp,1.0_dp,cs(1)))/real(n,dp)
    do i=2,n-1
      surv(i)=surv(i-1)*(1.0_dp-(merge(0.0_dp,1.0_dp,cs(i)))/real(n-i+1,dp))
    end do
    do i=1,n-1
      ints(i)=(ts(i+1)-ts(i))*surv(i)
    end do
    mr(n)=0.0_dp
    do i=n-1,1,-1
      if(surv(i)>0.0_dp)then
        mr(i)=(sum(ints(i:n-1)))/surv(i)
      else
        mr(i)=0.0_dp
      end if
    end do
    do i=1,n
      if(.not.cs(i))mr(i)=0.0_dp
      mrl(ord(i))=mr(i)
    end do
  end subroutine km_mean_residual_life
end module locfit_diagnostics
