! SPDX-License-Identifier: GPL-2.0-only
module ks_conditional
  use ks_kinds, only: dp
  use ks_bandwidth, only: hns_matrix
  use mvtnorm_probabilities, only: pmvnorm
  use mvtnorm_types, only: probability_result, probability_control
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_negative_inf, ieee_positive_inf
  implicit none
  private
  public :: kcde_model, fit_kcde, kcde_eval
  type :: kcde_model
    integer :: n=0,d=0
    real(dp),allocatable :: x(:,:),H(:,:),w(:)
  end type
contains
  subroutine fit_kcde(x,model,H,weights)
    real(dp),intent(in)::x(:,:)
    type(kcde_model),intent(out)::model
    real(dp),intent(in),optional::H(:,:),weights(:)
    real(dp)::sw
    model%n=size(x,1);model%d=size(x,2)
    if(model%n<=0.or.model%d<=0) error stop 'fit_kcde: empty data'
    allocate(model%x(model%n,model%d),model%H(model%d,model%d),model%w(model%n));model%x=x
    if(present(H))then
      if(any(shape(H)/=[model%d,model%d]))error stop 'fit_kcde: H shape';model%H=H
    else
      call hns_matrix(x,model%H)
    end if
    model%w=1.0_dp
    if(present(weights))then;if(size(weights)/=model%n)error stop 'fit_kcde: weights';model%w=weights;end if
    sw=sum(model%w);if(sw<=0.0_dp)error stop 'fit_kcde: nonpositive weight sum';model%w=model%w/sw
  end subroutine

  subroutine kcde_eval(model,points,values,upper_tail,control,errors)
    type(kcde_model),intent(in)::model
    real(dp),intent(in)::points(:,:)
    real(dp),intent(out)::values(size(points,1))
    logical,intent(in),optional::upper_tail
    type(probability_control),intent(in),optional::control
    real(dp),intent(out),optional::errors(size(points,1))
    type(probability_result)::res
    real(dp)::lo(model%d),up(model%d),ninf,pinf,err
    logical::upper
    integer::i,j
    if(size(points,2)/=model%d)error stop 'kcde_eval: point shape'
    ninf=ieee_value(0.0_dp,ieee_negative_inf);pinf=ieee_value(0.0_dp,ieee_positive_inf)
    upper=.false.;if(present(upper_tail))upper=upper_tail
    values=0.0_dp;if(present(errors))errors=0.0_dp
    do i=1,size(points,1)
      err=0.0_dp
      do j=1,model%n
        if(upper)then;lo=points(i,:);up=pinf
        else;lo=ninf;up=points(i,:)
        end if
        if(present(control))then;res=pmvnorm(lo,up,model%x(j,:),model%H,control)
        else;res=pmvnorm(lo,up,model%x(j,:),model%H)
        end if
        values(i)=values(i)+model%w(j)*res%value;err=err+model%w(j)*res%error
      end do
      if(present(errors))errors(i)=err
    end do
  end subroutine
end module ks_conditional
