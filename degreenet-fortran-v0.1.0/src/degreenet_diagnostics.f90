! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_diagnostics
  use degreenet_kinds, only : dp
  use degreenet_models, only : model_pmf
  implicit none
  private
  type, public :: mands_result
    integer, allocatable :: k(:)
    real(dp), allocatable :: ecdf(:), cdf(:), pmf(:)
    real(dp) :: statistic = 0.0_dp
  end type mands_result
  public :: modified_anderson_darling, hellinger_penalty, concentration_index

contains
  subroutine modified_anderson_darling(model,par,x,cutoff,cutabove,res)
    integer,intent(in)::model,x(:),cutoff,cutabove
    real(dp),intent(in)::par(:)
    type(mands_result),intent(out)::res
    integer::maxk,nk,i,j,n,kk
    real(dp)::cum,term
    n=count(x>=cutoff.and.x<=cutabove)
    if(n<=0)then;allocate(res%k(0),res%ecdf(0),res%cdf(0),res%pmf(0));return;end if
    maxk=maxval(x,mask=x>=cutoff.and.x<=cutabove);nk=maxk-cutoff+1
    allocate(res%k(nk),res%ecdf(nk),res%cdf(nk),res%pmf(nk));cum=0.0_dp
    do i=1,nk
      kk=cutoff+i-1;res%k(i)=kk;res%pmf(i)=model_pmf(model,par,kk,cutoff)
      cum=cum+res%pmf(i);res%cdf(i)=cum
      j=count(x>=cutoff.and.x<=kk.and.x<=cutabove);res%ecdf(i)=real(j,dp)/real(n,dp)
    end do
    res%statistic=0.0_dp
    do i=1,nk
      if(res%cdf(i)>0.0_dp.and.res%cdf(i)<1.0_dp)then
        term=(res%ecdf(i)-res%cdf(i))**2*res%pmf(i)/(res%cdf(i)*(1.0_dp-res%cdf(i)))
        res%statistic=res%statistic+real(n,dp)*term
      end if
    end do
  end subroutine modified_anderson_darling

  real(dp) function hellinger_penalty(model,par,x,cutoff,cutabove) result(h)
    integer,intent(in)::model,x(:),cutoff,cutabove
    real(dp),intent(in)::par(:)
    integer::k,n,nk
    real(dp)::emp,p,sum_mod,miss
    n=count(x>=cutoff.and.x<=cutabove);h=0.0_dp;sum_mod=0.0_dp
    if(n<=0)return
    do k=cutoff,maxval(x)
      nk=count(x==k);emp=real(nk,dp)/real(size(x),dp);p=model_pmf(model,par,k,cutoff)
      if(emp>0.0_dp)h=h+(sqrt(emp)-sqrt(max(p,0.0_dp)))**2
      if(emp>0.0_dp)sum_mod=sum_mod+p
    end do
    miss=max(0.0_dp,1.0_dp-sum_mod);h=h+0.5_dp*miss
  end function hellinger_penalty

  real(dp) function concentration_index(model,par,x,cutoff,maxk) result(c)
    integer,intent(in)::model,x(:),cutoff
    real(dp),intent(in)::par(:)
    integer,intent(in),optional::maxk
    integer::mk,k,n
    real(dp)::p0,p1,p2,m1,m2,pk
    mk=10000;if(present(maxk))mk=maxk;n=size(x)
    p0=real(count(x<cutoff),dp)/real(max(n,1),dp)
    p1=0.0_dp;p2=0.0_dp
    do k=1,cutoff-1
      p1=p1+real(k*count(x==k),dp)/real(max(n,1),dp)
      p2=p2+real(k*k*count(x==k),dp)/real(max(n,1),dp)
    end do
    m1=0.0_dp;m2=0.0_dp
    do k=cutoff,mk
      pk=model_pmf(model,par,k,cutoff);m1=m1+real(k,dp)*pk;m2=m2+real(k*k,dp)*pk
      if(k>cutoff+100.and.pk<1e-14_dp)exit
    end do
    c=(p1+m1*(1.0_dp-p0))/max(tiny(1.0_dp),p2+m2*(1.0_dp-p0)-p1-m1*(1.0_dp-p0))
  end function concentration_index
end module degreenet_diagnostics
