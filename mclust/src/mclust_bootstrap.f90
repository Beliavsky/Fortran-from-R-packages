! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_bootstrap
  use mclust_kinds, only : dp
  use mclust_types, only : mclust_fit, em_control
  use mclust_models, only : fit_model
  use mclust_simulation, only : simulate_fit
  implicit none
  private
  type, public :: bootstrap_lrt_result
    integer :: g0=0
    integer :: nboot=0
    real(dp) :: observed=0.0_dp
    real(dp) :: p_value=1.0_dp
    real(dp),allocatable :: statistic(:)
  end type
  type, public :: parameter_bootstrap_result
    integer :: nboot=0
    integer :: nonfit=0
    real(dp),allocatable :: pro(:,:)
    real(dp),allocatable :: mean(:,:,:)
    real(dp),allocatable :: variance(:,:,:,:)
  end type
  public :: bootstrap_lrt, mclust_parameter_bootstrap
contains
  subroutine bootstrap_lrt(x,model,g0,nboot,out,control,status)
    real(dp),intent(in)::x(:,:)
    character(len=*),intent(in)::model
    integer,intent(in)::g0,nboot
    type(bootstrap_lrt_result),intent(out)::out
    type(em_control),intent(in),optional::control
    integer,intent(out),optional::status
    type(mclust_fit)::f0,f1,b0,b1
    real(dp),allocatable::xb(:,:)
    integer,allocatable::comp(:)
    integer::b,info,attempts,max_attempts
    call fit_model(x,g0,model,f0,control); call fit_model(x,g0+1,model,f1,control)
    if(f0%status/=0 .or. f1%status/=0) then; if(present(status))status=-1; return; end if
    out%g0=g0; out%nboot=nboot; out%observed=2.0_dp*(f1%loglik-f0%loglik); allocate(out%statistic(nboot)); out%statistic=0.0_dp
    b=0; attempts=0; max_attempts=max(100,100*nboot)
    do while(b<nboot .and. attempts<max_attempts)
      attempts=attempts+1
      call simulate_fit(f0,size(x,1),xb,comp,info); if(info/=0) cycle
      call fit_model(xb,g0,model,b0,control); call fit_model(xb,g0+1,model,b1,control)
      if(b0%status/=0 .or. b1%status/=0) cycle
      b=b+1; out%statistic(b)=2.0_dp*(b1%loglik-b0%loglik)
    end do
    if(b<nboot) then
      if(present(status))status=1
      return
    end if
    out%p_value=(1.0_dp+real(count(out%statistic>=out%observed),dp))/real(nboot+1,dp)
    if(present(status))status=0
  end subroutine bootstrap_lrt

  subroutine mclust_parameter_bootstrap(fit,x,nboot,out,parametric,control,max_nonfit,status)
    type(mclust_fit),intent(in)::fit
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::nboot
    type(parameter_bootstrap_result),intent(out)::out
    logical,intent(in),optional::parametric
    type(em_control),intent(in),optional::control
    integer,intent(in),optional::max_nonfit
    integer,intent(out),optional::status
    logical::pb
    type(mclust_fit)::fb
    real(dp),allocatable::xb(:,:)
    integer,allocatable::comp(:)
    integer::b,i,j,info,nfail,mfail,n
    real(dp)::u
    pb=.false.; if(present(parametric))pb=parametric
    n=size(x,1); mfail=10*nboot; if(present(max_nonfit))mfail=max_nonfit
    out%nboot=nboot; allocate(out%pro(nboot,fit%g),out%mean(nboot,fit%d,fit%g),out%variance(nboot,fit%d,fit%d,fit%g))
    b=0; nfail=0
    do while(b<nboot .and. nfail<mfail)
      if(pb) then
        call simulate_fit(fit,n,xb,comp,info); if(info/=0) then; nfail=nfail+1; cycle; end if
      else
        allocate(xb(n,fit%d)); do i=1,n; call random_number(u); j=1+int(u*real(n,dp)); j=min(n,max(1,j)); xb(i,:)=x(j,:); end do
      end if
      call fit_model(xb,fit%g,fit%model_name,fb,control)
      if(fb%status/=0) then; nfail=nfail+1; deallocate(xb); if(allocated(comp))deallocate(comp); cycle; end if
      b=b+1; out%pro(b,:)=fb%pro; out%mean(b,:,:)=fb%mean; out%variance(b,:,:,:)=fb%sigma
      deallocate(xb); if(allocated(comp))deallocate(comp)
    end do
    out%nonfit=nfail
    if(b<nboot) then; if(present(status))status=1; else; if(present(status))status=0; end if
  end subroutine mclust_parameter_bootstrap
end module mclust_bootstrap
