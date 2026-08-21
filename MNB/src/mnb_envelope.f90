! SPDX-License-Identifier: GPL-2.0-or-later
module mnb_envelope
  use mnb_kinds, only : dp
  use mnb_types, only : mnb_fit_result,mnb_residual_result,mnb_envelope_result
  use mnb_core, only : fit_mnb
  use mnb_residuals, only : residuals_mnb,randomized_quantile_residuals
  use mnb_simulation, only : simulate_mnb
  implicit none
  private
  public :: envelope_mnb
contains
  function envelope_mnb(start,y,x,n,mi,residual_type,nsim,offset,source_compatible_quantile) result(out)
    real(dp),intent(in)::start(:),y(:),x(:,:)
    integer,intent(in)::n,mi,residual_type,nsim
    real(dp),intent(in),optional::offset(:)
    logical,intent(in),optional::source_compatible_quantile
    type(mnb_envelope_result)::out
    type(mnb_fit_result)::fit,fb
    type(mnb_residual_result)::r,rb
    real(dp),allocatable::tp(:),tb(:),e(:,:),rq(:)
    integer,allocatable::yb(:)
    integer::m,k
    logical::compat
    compat=.true.;if(present(source_compatible_quantile))compat=source_compatible_quantile
    fit=fit_mnb(start,y,x,n,mi,offset);r=residuals_mnb(start,y,x,n,mi,offset)
    call select_residual(residual_type,fit%par,y,r,tp)
    m=size(tp);allocate(e(m,nsim),yb(n*mi))
    do k=1,nsim
      call simulate_mnb(n,mi,x,fit%par,yb,offset)
      fb=fit_mnb(start,real(yb,dp),x,n,mi,offset);rb=residuals_mnb(start,real(yb,dp),x,n,mi,offset)
      if(residual_type==6 .and. compat)then
        allocate(rq(n));call randomized_quantile_residuals(fb%par,y,x,n,mi,rq,offset);allocate(tb(n));tb=rq;deallocate(rq)
      else
        call select_residual(residual_type,fb%par,real(yb,dp),rb,tb)
      end if
      call sort_ascending(tb);e(:,k)=tb;deallocate(tb)
    end do
    allocate(out%lower(m),out%mean(m),out%upper(m),out%residual(m));out%residual=tp
    do k=1,m;out%lower(k)=minval(e(k,:));out%upper(k)=maxval(e(k,:));out%mean(k)=sum(e(k,:))/real(nsim,dp);end do
  contains
    subroutine select_residual(kind,par,yy,rr,v)
      integer,intent(in)::kind;real(dp),intent(in)::par(:),yy(:);type(mnb_residual_result),intent(in)::rr
      real(dp),allocatable,intent(out)::v(:)
      select case(kind)
      case(1);allocate(v(size(rr%weighted)));v=rr%weighted
      case(2);allocate(v(size(rr%standardized_weighted)));v=rr%standardized_weighted
      case(3);allocate(v(size(rr%pearson)));v=rr%pearson
      case(4);allocate(v(size(rr%standardized_pearson)));v=rr%standardized_pearson
      case(5);allocate(v(size(rr%deviance)));v=rr%deviance
      case(6);allocate(v(n));call randomized_quantile_residuals(par,yy,x,n,mi,v,offset)
      case default;error stop 'envelope_mnb: residual type must be 1..6'
      end select
    end subroutine
    subroutine sort_ascending(v)
      real(dp),intent(inout)::v(:);real(dp)::key;integer::i,j
      do i=2,size(v);key=v(i);j=i-1;do while(j>=1);if(v(j)<=key)exit;v(j+1)=v(j);j=j-1;end do;v(j+1)=key;end do
    end subroutine
  end function envelope_mnb
end module mnb_envelope
