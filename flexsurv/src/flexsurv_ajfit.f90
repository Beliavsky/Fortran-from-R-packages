! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_ajfit
  use flexsurv_kinds, only : dp
  use flexsurv_multistate, only : flexsurv_transition, pmatrix_flexsurv
  implicit none
  private

  type, public :: ajfit_result
    real(dp), allocatable :: time(:)
    real(dp), allocatable :: probability(:,:)
  end type ajfit_result

  type, public :: ajfit_comparison
    real(dp), allocatable :: time(:)
    real(dp), allocatable :: aj(:,:)
    real(dp), allocatable :: parametric(:,:)
  end type ajfit_comparison

  public :: ajfit_transition_data, ajfit_compare_parametric

contains

  function ajfit_transition_data(time,status,from,to,nstate,start_state) result(res)
    real(dp),intent(in)::time(:)
    integer,intent(in)::status(:),from(:),to(:),nstate,start_state
    type(ajfit_result)::res
    real(dp),allocatable::et(:),q(:,:),p(:),pnew(:)
    real(dp)::risk,d
    integer::ne,k,r,a,b,j

    call unique_event_times(time,status,et);ne=size(et)
    allocate(res%time(ne+1),res%probability(nstate,ne+1),p(nstate),pnew(nstate),q(nstate,nstate))
    res%time(1)=0.0_dp;res%probability(:,1)=0.0_dp;p=0.0_dp;p(start_state)=1.0_dp
    res%probability(:,1)=p
    do k=1,ne
      q=0.0_dp
      do r=1,size(time)
        if(status(r)/=1)cycle
        if(.not.same_time(time(r),et(k)))cycle
        a=from(r);b=to(r)
        if(a<1.or.a>nstate.or.b<1.or.b>nstate)cycle
        risk=0.0_dp
        do j=1,size(time)
          if(from(j)==a.and.to(j)==b.and.time(j)>=et(k))risk=risk+1.0_dp
        end do
        if(risk<=0.0_dp)cycle
        d=1.0_dp
        q(a,b)=q(a,b)+d/risk
      end do
      do a=1,nstate
        q(a,a)=-sum(q(a,:))+q(a,a)
      end do
      pnew=p+matmul(p,q);p=pnew
      where(abs(p)<100.0_dp*epsilon(1.0_dp))p=0.0_dp
      res%time(k+1)=et(k);res%probability(:,k+1)=p
    end do
  end function ajfit_transition_data

  function ajfit_compare_parametric(trans,time,status,from,to,nstate,start_state) result(out)
    type(flexsurv_transition),intent(in)::trans(:)
    real(dp),intent(in)::time(:)
    integer,intent(in)::status(:),from(:),to(:),nstate,start_state
    type(ajfit_comparison)::out
    type(ajfit_result)::aj
    real(dp),allocatable::pmat(:,:,:)
    integer::st
    aj=ajfit_transition_data(time,status,from,to,nstate,start_state)
    call pmatrix_flexsurv(trans,nstate,aj%time,pmat,st)
    allocate(out%time(size(aj%time)),out%aj(nstate,size(aj%time)))
    out%time=aj%time;out%aj=aj%probability
    allocate(out%parametric(nstate,size(aj%time)));out%parametric=pmat(start_state,:,:)
  end function ajfit_compare_parametric

  subroutine unique_event_times(time,status,u)
    real(dp),intent(in)::time(:)
    integer,intent(in)::status(:)
    real(dp),allocatable,intent(out)::u(:)
    real(dp),allocatable::a(:),tmp(:)
    real(dp)::key
    integer::n,i,j,m
    n=count(status==1);allocate(a(n));m=0
    do i=1,size(time)
      if(status(i)==1)then;m=m+1;a(m)=time(i);end if
    end do
    do i=2,n
      key=a(i);j=i-1
      do while(j>=1)
        if(a(j)<=key)exit
        a(j+1)=a(j);j=j-1
      end do
      a(j+1)=key
    end do
    allocate(tmp(n));m=0
    do i=1,n
      if(m==0)then
        m=1;tmp(m)=a(i)
      else if(.not.same_time(a(i),tmp(m)))then
        m=m+1;tmp(m)=a(i)
      end if
    end do
    allocate(u(m));if(m>0)u=tmp(1:m)
  end subroutine unique_event_times

  logical function same_time(a,b) result(eq)
    real(dp),intent(in)::a,b
    eq=abs(a-b)<=32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(a),abs(b))
  end function same_time

end module flexsurv_ajfit
