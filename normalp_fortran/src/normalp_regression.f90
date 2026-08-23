module normalp_regression
  use normalp_special, only: dp
  use normalp_estimation, only: estimatep
  implicit none
  private
  public :: lmp_fit, lmp_result
  type :: lmp_result
    real(dp), allocatable :: coef(:), fitted(:), residuals(:)
    real(dp) :: p=2.0_dp, sigma=0.0_dp
    integer :: iterations=0, no_conv=0
  end type
contains
  subroutine lmp_fit(x,y,res,p_fixed,tol,max_iter)
    real(dp),intent(in)::x(:,:),y(:); type(lmp_result),intent(out)::res
    real(dp),intent(in),optional::p_fixed,tol; integer,intent(in),optional::max_iter
    integer::n,k,it,mx; real(dp)::p,oldp,tolv; real(dp),allocatable::b(:),oldb(:),r(:)
    n=size(y); k=size(x,2); mx=100; if(present(max_iter)) mx=max_iter; tolv=1.0e-4_dp; if(present(tol)) tolv=tol
    allocate(b(k),oldb(k),r(n)); call ols_start(x,y,b); r=y-matmul(x,b)
    if(present(p_fixed)) then
      p=p_fixed; call lp_reg_opt(x,y,p,b); it=1
    else
      p=estimatep(r,sum(r)/real(n,dp),2.0_dp)
      do it=1,mx
        oldp=p; oldb=b; call lp_reg_opt(x,y,p,b); r=y-matmul(x,b); p=estimatep(r,sum(r)/real(n,dp),p)
        if(abs(p-oldp)<=tolv .and. sum(abs(b-oldb))<=tolv) exit
      end do
      if(it>mx) then; res%no_conv=1; it=mx; end if
    end if
    allocate(res%coef(k),res%fitted(n),res%residuals(n)); res%coef=b; res%fitted=matmul(x,b); res%residuals=y-res%fitted; res%p=p
    res%sigma=(sum(abs(res%residuals)**p)/real(max(1,n-k-merge(1,0,.not.present(p_fixed))),dp))**(1.0_dp/p)
    res%iterations=it
  end subroutine

  subroutine ols_start(x,y,b)
    real(dp),intent(in)::x(:,:),y(:); real(dp),intent(out)::b(:)
    real(dp),allocatable::a(:,:),rhs(:); integer::k
    k=size(x,2); allocate(a(k,k),rhs(k)); a=matmul(transpose(x),x); rhs=matmul(transpose(x),y); call solve_linear(a,rhs,b)
  end subroutine

  subroutine lp_reg_opt(x,y,p,b)
    real(dp),intent(in)::x(:,:),y(:),p; real(dp),intent(inout)::b(:)
    real(dp)::step,best,val; real(dp),allocatable::bt(:); integer::i,it
    allocate(bt(size(b))); best=obj(b); step=max(0.1_dp,maxval(abs(b))*0.25_dp)
    do it=1,400
      do i=1,size(b)
        bt=b; bt(i)=bt(i)+step; val=obj(bt); if(val<best) then; b=bt; best=val; cycle; end if
        bt=b; bt(i)=bt(i)-step; val=obj(bt); if(val<best) then; b=bt; best=val; end if
      end do
      step=step*0.92_dp; if(step<1.0e-8_dp) exit
    end do
  contains
    function obj(bb) result(v)
      real(dp),intent(in)::bb(:); real(dp)::v
      v=sum(abs(y-matmul(x,bb))**p)
    end function
  end subroutine

  subroutine solve_linear(a,b,x)
    real(dp),intent(in)::a(:,:),b(:); real(dp),intent(out)::x(:)
    real(dp),allocatable::m(:,:),rhs(:); real(dp)::fac,tmp; integer::n,i,j,k,piv
    n=size(b); m=a; rhs=b
    do k=1,n-1
      piv=k; do i=k+1,n; if(abs(m(i,k))>abs(m(piv,k))) piv=i; end do
      if(piv/=k) then; do j=k,n; tmp=m(k,j);m(k,j)=m(piv,j);m(piv,j)=tmp;end do; tmp=rhs(k);rhs(k)=rhs(piv);rhs(piv)=tmp; end if
      do i=k+1,n; fac=m(i,k)/m(k,k); m(i,k:n)=m(i,k:n)-fac*m(k,k:n); rhs(i)=rhs(i)-fac*rhs(k); end do
    end do
    x(n)=rhs(n)/m(n,n); do i=n-1,1,-1; x(i)=(rhs(i)-sum(m(i,i+1:n)*x(i+1:n)))/m(i,i); end do
  end subroutine
end module normalp_regression
