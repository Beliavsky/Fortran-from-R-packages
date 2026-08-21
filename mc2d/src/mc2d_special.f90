! SPDX-License-Identifier: GPL-2.0-or-later
module mc2d_special
  use mc2d_kinds, only : dp, nan_dp
  use mvtnorm_special, only : regularized_beta, normal_cdf, normal_quantile
  implicit none
  private
  public :: beta_pdf, beta_cdf, beta_quantile
  public :: ncbeta_pdf, ncbeta_cdf, ncbeta_quantile
  public :: normal_cdf, normal_quantile
contains
  real(dp) function beta_pdf(x,a,b,log_density) result(v)
    real(dp),intent(in)::x,a,b
    logical,intent(in),optional::log_density
    logical::llog
    real(dp)::lv
    llog=.false.; if(present(log_density)) llog=log_density
    if(a<=0.0_dp .or. b<=0.0_dp) then; v=nan_dp(); return; end if
    if(x<0.0_dp .or. x>1.0_dp) then
      if(llog) v=-huge(1.0_dp); if(.not.llog) v=0.0_dp; return
    end if
    if(x==0.0_dp) then
      if(a<1.0_dp) then; lv=huge(1.0_dp); else if(a==1.0_dp) then; lv=log(b); else; lv=-huge(1.0_dp); end if
    else if(x==1.0_dp) then
      if(b<1.0_dp) then; lv=huge(1.0_dp); else if(b==1.0_dp) then; lv=log(a); else; lv=-huge(1.0_dp); end if
    else
      lv=(a-1.0_dp)*log(x)+(b-1.0_dp)*log(1.0_dp-x)-log_gamma(a)-log_gamma(b)+log_gamma(a+b)
    end if
    if(llog) then; v=lv; else; v=exp(lv); end if
  end function beta_pdf

  real(dp) function beta_cdf(x,a,b) result(v)
    real(dp),intent(in)::x,a,b
    if(a<=0.0_dp .or. b<=0.0_dp) then; v=nan_dp(); else; v=regularized_beta(x,a,b); end if
  end function beta_cdf

  real(dp) function beta_quantile(p,a,b) result(x)
    real(dp),intent(in)::p,a,b
    real(dp)::lo,hi,mid,f,pdf
    integer::i
    if(p<0.0_dp .or. p>1.0_dp .or. a<=0.0_dp .or. b<=0.0_dp) then; x=nan_dp(); return; end if
    if(p==0.0_dp) then; x=0.0_dp; return; end if
    if(p==1.0_dp) then; x=1.0_dp; return; end if
    lo=0.0_dp; hi=1.0_dp; x=a/(a+b)
    do i=1,100
      f=beta_cdf(x,a,b)-p
      if(f>0.0_dp) hi=x; if(f<0.0_dp) lo=x
      pdf=beta_pdf(x,a,b)
      if(pdf>tiny(1.0_dp) .and. pdf<huge(1.0_dp)/10.0_dp) then; mid=x-f/pdf; else; mid=0.5_dp*(lo+hi); end if
      if(mid<=lo .or. mid>=hi) mid=0.5_dp*(lo+hi)
      if(abs(mid-x)<=16*epsilon(1.0_dp)*max(1.0_dp,abs(x))) then; x=mid; exit; end if
      x=mid
    end do
  end function beta_quantile

  real(dp) function ncbeta_cdf(x,a,b,ncp) result(v)
    real(dp),intent(in)::x,a,b,ncp
    real(dp)::lam,w,sumw,term
    integer::k
    if(ncp<0.0_dp) then; v=nan_dp(); return; end if
    if(ncp==0.0_dp) then; v=beta_cdf(x,a,b); return; end if
    lam=0.5_dp*ncp; w=exp(-lam); sumw=0.0_dp; v=0.0_dp; k=0
    do
      term=w*beta_cdf(x,a+real(k,dp),b); v=v+term; sumw=sumw+w
      if(1.0_dp-sumw < 16*epsilon(1.0_dp) .and. k>lam+10*sqrt(lam+1.0_dp)) exit
      k=k+1; w=w*lam/real(k,dp)
      if(k>100000) exit
    end do
    v=min(1.0_dp,max(0.0_dp,v))
  end function ncbeta_cdf

  real(dp) function ncbeta_pdf(x,a,b,ncp,log_density) result(v)
    real(dp),intent(in)::x,a,b,ncp
    logical,intent(in),optional::log_density
    real(dp)::lam,w,sumw,d
    integer::k
    logical::llog
    llog=.false.; if(present(log_density)) llog=log_density
    if(ncp<0.0_dp) then; v=nan_dp(); return; end if
    if(ncp==0.0_dp) then; v=beta_pdf(x,a,b,llog); return; end if
    lam=0.5_dp*ncp; w=exp(-lam); sumw=0.0_dp; d=0.0_dp; k=0
    do
      d=d+w*beta_pdf(x,a+real(k,dp),b); sumw=sumw+w
      if(1.0_dp-sumw < 16*epsilon(1.0_dp) .and. k>lam+10*sqrt(lam+1.0_dp)) exit
      k=k+1; w=w*lam/real(k,dp); if(k>100000) exit
    end do
    if(llog) then; if(d>0) v=log(d); if(d<=0) v=-huge(1.0_dp); else; v=d; end if
  end function ncbeta_pdf

  real(dp) function ncbeta_quantile(p,a,b,ncp) result(x)
    real(dp),intent(in)::p,a,b,ncp
    real(dp)::lo,hi,mid
    integer::i
    if(p<0.0_dp .or. p>1.0_dp .or. a<=0.0_dp .or. b<=0.0_dp .or. ncp<0.0_dp) then; x=nan_dp(); return; end if
    if(p==0.0_dp) then; x=0.0_dp; return; end if
    if(p==1.0_dp) then; x=1.0_dp; return; end if
    lo=0.0_dp; hi=1.0_dp
    do i=1,100
      mid=0.5_dp*(lo+hi)
      if(ncbeta_cdf(mid,a,b,ncp)<p) then; lo=mid; else; hi=mid; end if
    end do
    x=0.5_dp*(lo+hi)
  end function ncbeta_quantile
end module mc2d_special
