! SPDX-License-Identifier: GPL-2.0-or-later
module mnb_optimizer
  use mnb_kinds, only : dp
  use mnb_math, only : numerical_gradient
  implicit none
  private
  public :: bfgs_maximize
  abstract interface
    function objective_fn(x) result(f)
      import dp
      real(dp),intent(in)::x(:)
      real(dp)::f
    end function objective_fn
  end interface
contains
  subroutine bfgs_maximize(fn,x0,x,fval,iters,status,maxit,tol)
    procedure(objective_fn)::fn
    real(dp),intent(in)::x0(:)
    real(dp),intent(out)::x(:),fval
    integer,intent(out)::iters,status
    integer,intent(in),optional::maxit
    real(dp),intent(in),optional::tol
    integer::n,i,mx
    real(dp)::gtol,alpha,fnew,ys,rho
    real(dp),allocatable::h(:,:),g(:),gnew(:),p(:),s(:),y(:),xnew(:),eye(:,:)
    n=size(x0);mx=500;if(present(maxit))mx=maxit;gtol=1.0e-6_dp;if(present(tol))gtol=tol
    allocate(h(n,n),g(n),gnew(n),p(n),s(n),y(n),xnew(n),eye(n,n));h=0.0_dp;eye=0.0_dp
    do i=1,n;h(i,i)=1.0_dp;eye(i,i)=1.0_dp;end do
    x=x0; fval=fn(x); call grad_neg(x,g)
    status=1
    do iters=1,mx
      if(maxval(abs(g))<gtol)then;status=0;exit;end if
      p=-matmul(h,g); alpha=1.0_dp
      do
        xnew=x+alpha*p
        if(xnew(1)<=1.0e-10_dp)then;alpha=alpha*0.5_dp;cycle;end if
        fnew=fn(xnew)
        if(fnew> -huge(1.0_dp)/4.0_dp .and. -fnew <= -fval + 1.0e-4_dp*alpha*dot_product(g,p)) exit
        alpha=alpha*0.5_dp
        if(alpha<1.0e-10_dp)exit
      end do
      if(alpha<1.0e-10_dp)exit
      call grad_neg(xnew,gnew);s=xnew-x;y=gnew-g;ys=dot_product(y,s)
      if(ys>1.0e-14_dp)then
        rho=1.0_dp/ys
        h=matmul(matmul(eye-rho*outer(s,y),h),eye-rho*outer(y,s))+rho*outer(s,s)
      else
        h=eye
      end if
      x=xnew;g=gnew;fval=fnew
    end do
  contains
    subroutine grad_neg(xx,gg)
      real(dp),intent(in)::xx(:);real(dp),intent(out)::gg(:)
      call numerical_gradient(negfn,xx,gg)
    end subroutine
    real(dp) function negfn(xx) result(v)
      real(dp),intent(in)::xx(:);v=-fn(xx)
    end function
    pure function outer(a,b) result(c)
      real(dp),intent(in)::a(:),b(:);real(dp)::c(size(a),size(b));integer::ii
      do ii=1,size(a);c(ii,:)=a(ii)*b;end do
    end function
  end subroutine bfgs_maximize
end module mnb_optimizer
