! SPDX-License-Identifier: GPL-2.0-only
module optimx_example_functions
  use optimx_kinds,only:dp
  implicit none
  private
  public::rosenbrock_callback,quadratic_callback,quartic_callback
contains
  subroutine rosenbrock_callback(x,value,gradient,hessian,need_gradient,need_hessian,status)
    real(dp),intent(in)::x(:);real(dp),intent(out)::value
    real(dp),intent(inout)::gradient(:),hessian(:,:)
    logical,intent(in)::need_gradient,need_hessian;integer,intent(out)::status
    value=(1.0_dp-x(1))**2+100.0_dp*(x(2)-x(1)**2)**2;status=0
    if(need_gradient)then
      gradient=0.0_dp;gradient(1)=-2.0_dp*(1.0_dp-x(1))-400.0_dp*x(1)*(x(2)-x(1)**2)
      gradient(2)=200.0_dp*(x(2)-x(1)**2)
    end if
    if(need_hessian)then
      hessian=0.0_dp;hessian(1,1)=2.0_dp-400.0_dp*x(2)+1200.0_dp*x(1)**2
      hessian(1,2)=-400.0_dp*x(1);hessian(2,1)=hessian(1,2);hessian(2,2)=200.0_dp
    end if
  end subroutine rosenbrock_callback
  subroutine quadratic_callback(x,value,gradient,hessian,need_gradient,need_hessian,status)
    real(dp),intent(in)::x(:);real(dp),intent(out)::value
    real(dp),intent(inout)::gradient(:),hessian(:,:)
    logical,intent(in)::need_gradient,need_hessian;integer,intent(out)::status
    integer :: i
    value=sum((x-[(real(i,dp),i=1,size(x))])**2);status=0
    if(need_gradient)gradient=2.0_dp*(x-[(real(i,dp),i=1,size(x))])
    if(need_hessian)then;hessian=0.0_dp;block;integer::i;do i=1,size(x);hessian(i,i)=2.0_dp;end do;end block;end if
  end subroutine quadratic_callback
  subroutine quartic_callback(x,value,gradient,hessian,need_gradient,need_hessian,status)
    real(dp),intent(in)::x(:);real(dp),intent(out)::value
    real(dp),intent(inout)::gradient(:),hessian(:,:)
    logical,intent(in)::need_gradient,need_hessian;integer,intent(out)::status
    value=(x(1)*x(1)-1.0_dp)**2+0.1_dp*(x(2)-2.0_dp)**2;status=0
    if(need_gradient)then;gradient=0.0_dp;gradient(1)=4.0_dp*x(1)*(x(1)*x(1)-1.0_dp);gradient(2)=0.2_dp*(x(2)-2.0_dp);end if
    if(need_hessian)then;hessian=0.0_dp;hessian(1,1)=12.0_dp*x(1)*x(1)-4.0_dp;hessian(2,2)=0.2_dp;end if
  end subroutine quartic_callback
end module optimx_example_functions
