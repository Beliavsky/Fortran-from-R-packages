! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_simulation
  use degreenet_kinds, only : dp
  use degreenet_rng, only : runif01
  use degreenet_models, only : model_pmf
  use degreenet_fit, only : fit_result, fit_degree_model
  implicit none
  private
  public :: sample_model, bootstrap_degree_model, bootstrap_resample_fit

contains
  subroutine sample_model(model,par,n,sample,cutoff,maxdeg)
    integer,intent(in)::model,n,cutoff,maxdeg
    real(dp),intent(in)::par(:)
    integer,intent(out)::sample(n)
    real(dp),allocatable::cdf(:)
    real(dp)::s,u
    integer::k,i,j
    allocate(cdf(cutoff:maxdeg));s=0.0_dp
    do k=cutoff,maxdeg;s=s+model_pmf(model,par,k,cutoff);cdf(k)=s;end do
    if(s<=0.0_dp)then;sample=cutoff;return;end if
    cdf=cdf/s
    do i=1,n
      u=runif01();j=cutoff
      do while(j<maxdeg.and.u>cdf(j));j=j+1;end do
      sample(i)=j
    end do
  end subroutine sample_model

  subroutine bootstrap_degree_model(model,theta,n,b,cutoff,cutabove,maxdeg,start,estimates, &
      lower,upper)
    integer,intent(in)::model,n,b,cutoff,cutabove,maxdeg
    real(dp),intent(in)::theta(:),start(:)
    real(dp),intent(out)::estimates(size(start),b)
    real(dp),intent(in),optional::lower(:),upper(:)
    integer,allocatable::x(:)
    type(fit_result)::fit
    integer::j
    allocate(x(n))
    do j=1,b
      call sample_model(model,theta,n,x,cutoff,maxdeg)
      call fit_degree_model(model,x,cutoff,cutabove,start,fit,lower,upper,maxit=1500)
      estimates(:,j)=fit%theta
    end do
  end subroutine bootstrap_degree_model

  subroutine bootstrap_resample_fit(model,x,b,cutoff,cutabove,start,estimates,lower,upper)
    integer,intent(in)::model,x(:),b,cutoff,cutabove
    real(dp),intent(in)::start(:)
    real(dp),intent(out)::estimates(size(start),b)
    real(dp),intent(in),optional::lower(:),upper(:)
    integer,allocatable::xs(:)
    integer::i,j,ix
    type(fit_result)::fit
    allocate(xs(size(x)))
    do j=1,b
      do i=1,size(x)
        ix=1+int(runif01()*real(size(x),dp))
        ix=min(size(x),max(1,ix));xs(i)=x(ix)
      end do
      call fit_degree_model(model,xs,cutoff,cutabove,start,fit,lower,upper,maxit=1500)
      estimates(:,j)=fit%theta
    end do
  end subroutine bootstrap_resample_fit
end module degreenet_simulation
