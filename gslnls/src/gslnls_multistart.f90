! SPDX-License-Identifier: LGPL-3.0-only
module gslnls_multistart
  use gslnls_kinds, only : dp
  use gslnls_types
  use gslnls_core, only : fit_nls
  use gslnls_linalg, only : determinant_jtj
  implicit none
  private
  public :: fit_nls_multistart, halton_point

contains

  subroutine fit_nls_multistart(model, y, ranges, result, control, jac, fvv, lower, upper, weights, loss)
    procedure(nls_model) :: model
    real(dp), intent(in) :: y(:), ranges(:,:)
    type(multistart_result), intent(out) :: result
    type(nls_control), intent(in), optional :: control
    procedure(nls_jacobian), optional :: jac
    procedure(nls_fvv), optional :: fvv
    real(dp), intent(in), optional :: lower(:), upper(:), weights(:)
    type(nls_loss), intent(in), optional :: loss
    type(nls_control) :: ctl, cheap_ctl, full_ctl
    type(nls_result) :: tmp
    real(dp), allocatable :: points(:,:), ssr(:), best(:)
    integer, allocatable :: keep(:), keep2(:), order(:)
    real(dp) :: best_ssr, detj
    integer :: p, npt, major, i, k, idx, nsp, nwsp, q, sample_index
    logical :: improved

    ctl = nls_control(); if (present(control)) ctl = control
    p = size(ranges,2); npt = ctl%mstart_n
    if (size(ranges,1) /= 2 .or. any(ranges(1,:) > ranges(2,:)) .or. p < 1 .or. npt < 1) then
      result%fit%status = NLS_BAD_INPUT
      return
    end if
    q = min(max(1,ctl%mstart_q),npt)
    allocate(points(npt,p),ssr(npt),best(p),keep(npt),keep2(npt),order(npt))
    points = 0.0_dp; ssr = huge(1.0_dp); keep = 0
    best = 0.5_dp*(ranges(1,:)+ranges(2,:)); best_ssr=huge(1.0_dp)
    nsp=0; nwsp=0; sample_index=0
    cheap_ctl=ctl; cheap_ctl%maxiter=ctl%mstart_p; cheap_ctl%gtol=1.0e-3_dp; cheap_ctl%store_trace=.false.
    full_ctl=ctl; full_ctl%store_trace=.false.

    do major=1,ctl%mstart_maxstart
      do i=1,npt
        if(keep(i)==0) then
          sample_index=sample_index+1
          call sample_range(sample_index,ranges,points(i,:))
        end if
        call run_fit_dispatch(model,y,points(i,:),cheap_ctl,tmp,jac,fvv,lower,upper,weights,loss)
        result%starts_evaluated=result%starts_evaluated+1
        if(tmp%status==NLS_SUCCESS .or. tmp%status==NLS_MAXITER) then
          points(i,:)=tmp%par; ssr(i)=tmp%ssr
          if(allocated(tmp%jacobian)) then
            detj=determinant_jtj(tmp%jacobian)
            if(detj<=1.0e-12_dp .and. tmp%ssr>2.0_dp*ctl%ftol) ssr(i)=huge(1.0_dp)
          end if
        else
          ssr(i)=huge(1.0_dp)
        end if
      end do
      call order_real(ssr,order)
      keep2=0
      do k=1,q
        idx=order(k)
        if(ssr(idx)<huge(1.0_dp)/4.0_dp) keep2(idx)=keep(idx)+1
      end do
      keep=keep2
      improved=.false.
      do k=1,q
        idx=order(k)
        if(keep(idx)>=ctl%mstart_s .and. ssr(idx)<huge(1.0_dp)/4.0_dp) then
          keep(idx)=0; result%local_searches=result%local_searches+1
          call run_fit_dispatch(model,y,points(idx,:),full_ctl,tmp,jac,fvv,lower,upper,weights,loss)
          if(tmp%status==NLS_SUCCESS .or. tmp%status==NLS_MAXITER) then
            if(tmp%ssr < 0.99_dp*best_ssr) then
              best_ssr=tmp%ssr; best=tmp%par; result%fit=tmp
              nsp=nsp+1; nwsp=0; improved=.true.
            else
              nwsp=nwsp+1
            end if
          else
            nwsp=nwsp+1
          end if
        end if
      end do
      if(.not.improved .and. nsp>0) nwsp=nwsp+1
      result%major_iterations=major
      if(nsp>=ctl%mstart_minsp .and. real(nwsp,dp) > ctl%mstart_r + sqrt(ctl%mstart_r)*real(nsp,dp)) exit
    end do

    if(best_ssr>=huge(1.0_dp)/4.0_dp) then
      call run_fit_dispatch(model,y,best,full_ctl,result%fit,jac,fvv,lower,upper,weights,loss)
    end if
    result%best_start=best

  end subroutine fit_nls_multistart


  subroutine run_fit_dispatch(model,y,x0,cc,rr,jac,fvv,lower,upper,weights,loss)
    procedure(nls_model) :: model
    real(dp), intent(in) :: y(:), x0(:)
    type(nls_control), intent(in) :: cc
    type(nls_result), intent(out) :: rr
    procedure(nls_jacobian), optional :: jac
    procedure(nls_fvv), optional :: fvv
    real(dp), intent(in), optional :: lower(:), upper(:), weights(:)
    type(nls_loss), intent(in), optional :: loss

    if (present(jac)) then
      if (present(fvv)) then
        if (present(loss)) then
          call fit_nls(model,y,x0,rr,cc,jac=jac,fvv=fvv,lower=lower,upper=upper,weights=weights,loss=loss)
        else
          call fit_nls(model,y,x0,rr,cc,jac=jac,fvv=fvv,lower=lower,upper=upper,weights=weights)
        end if
      else
        if (present(loss)) then
          call fit_nls(model,y,x0,rr,cc,jac=jac,lower=lower,upper=upper,weights=weights,loss=loss)
        else
          call fit_nls(model,y,x0,rr,cc,jac=jac,lower=lower,upper=upper,weights=weights)
        end if
      end if
    else
      if (present(loss)) then
        call fit_nls(model,y,x0,rr,cc,lower=lower,upper=upper,weights=weights,loss=loss)
      else
        call fit_nls(model,y,x0,rr,cc,lower=lower,upper=upper,weights=weights)
      end if
    end if
  end subroutine run_fit_dispatch

  subroutine sample_range(index,ranges,x)
    integer,intent(in)::index
    real(dp),intent(in)::ranges(:,:)
    real(dp),intent(out)::x(:)
    real(dp),allocatable::u(:)
    allocate(u(size(x))); call halton_point(index,size(x),u)
    x=ranges(1,:)+(ranges(2,:)-ranges(1,:))*u
  end subroutine sample_range

  subroutine halton_point(index,dim,u)
    integer,intent(in)::index,dim
    real(dp),intent(out)::u(dim)
    integer::j
    do j=1,dim
      u(j)=radical_inverse(index,nth_prime(j))
    end do
  end subroutine halton_point

  pure real(dp) function radical_inverse(index,base) result(v)
    integer,intent(in)::index,base
    integer::n
    real(dp)::f
    n=index; v=0.0_dp; f=1.0_dp/real(base,dp)
    do while(n>0)
      v=v+f*real(mod(n,base),dp); n=n/base; f=f/real(base,dp)
    end do
  end function radical_inverse

  integer function nth_prime(k) result(p)
    integer,intent(in)::k
    integer::count,cand,d
    logical::prime
    count=0; cand=1
    do while(count<k)
      cand=cand+1; prime=.true.
      do d=2,int(sqrt(real(cand,dp)))
        if(mod(cand,d)==0) then; prime=.false.; exit; end if
      end do
      if(prime) count=count+1
    end do
    p=cand
  end function nth_prime

  subroutine order_real(x,ord)
    real(dp),intent(in)::x(:)
    integer,intent(out)::ord(:)
    integer::i,j,t
    do i=1,size(x); ord(i)=i; end do
    do i=2,size(x)
      t=ord(i); j=i-1
      do while(j>=1)
        if(x(ord(j))<=x(t)) exit
        ord(j+1)=ord(j); j=j-1
      end do
      ord(j+1)=t
    end do
  end subroutine order_real

end module gslnls_multistart
