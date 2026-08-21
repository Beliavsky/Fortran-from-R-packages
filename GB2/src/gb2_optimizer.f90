! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_optimizer
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use gb2_kinds, only : dp
  implicit none
  private
  public :: optimization_result, bfgs_minimize, invert_matrix
  type :: optimization_result
    real(dp), allocatable :: par(:)
    real(dp) :: value=huge(1.0_dp)
    integer :: iterations=0
    logical :: converged=.false.
  end type optimization_result
  abstract interface
    function objective_fun(x) result(v)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: v
    end function objective_fun
    subroutine gradient_fun(x,g)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
    end subroutine gradient_fun
  end interface
contains
  subroutine bfgs_minimize(fn,gr,x0,result,maxiter,tol)
    procedure(objective_fun) :: fn
    procedure(gradient_fun) :: gr
    real(dp), intent(in) :: x0(:)
    type(optimization_result), intent(out) :: result
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol
    integer :: n,it,nmax,i
    real(dp) :: eps,f,fnw,alpha,rho,ys
    real(dp), allocatable :: x(:),xn(:),g(:),gn(:),p(:),s(:),y(:),h(:,:),eye(:,:),tmp(:,:)
    n=size(x0)
    nmax=500
    if(present(maxiter)) nmax=max(1,maxiter)
    eps=1.0e-8_dp
    if(present(tol)) eps=max(tol,10.0_dp*epsilon(1.0_dp))
    allocate(x(n),xn(n),g(n),gn(n),p(n),s(n),y(n),h(n,n),eye(n,n),tmp(n,n))
    eye=0.0_dp
    do i=1,n
    eye(i,i)=1.0_dp
    end do
    h=eye
    x=x0
    f=fn(x)
    call gr(x,g)
    result%converged=.false.
    it=0
    do it=1,nmax
      if(maxval(abs(g))<=eps*(1.0_dp+abs(f))) then
      result%converged=.true.
      exit
      end if
      p=-matmul(h,g)
      if(dot_product(p,g)>=0.0_dp) then
      p=-g
      h=eye
      end if
      alpha=1.0_dp
      do
        xn=x+alpha*p
        fnw=fn(xn)
        if(ieee_is_finite(fnw) .and. fnw<=f+1.0e-4_dp*alpha*dot_product(g,p)) exit
        alpha=0.5_dp*alpha
        if(alpha<1.0e-12_dp) exit
      end do
      if(alpha<1.0e-12_dp .or. .not.ieee_is_finite(fnw)) exit
      call gr(xn,gn)
      s=xn-x
      y=gn-g
      ys=dot_product(y,s)
      if(ys>1.0e-14_dp*sqrt(dot_product(y,y)*dot_product(s,s))) then
        rho=1.0_dp/ys
        tmp=eye-rho*outer(s,y)
        h=matmul(tmp,matmul(h,transpose(tmp)))+rho*outer(s,s)
      else
        h=eye
      end if
      x=xn
      g=gn
      f=fnw
      if(sqrt(dot_product(s,s))<=eps*(1.0_dp+sqrt(dot_product(x,x)))) then
        result%converged=.true.
        exit
      end if
    end do
    allocate(result%par(n))
    result%par=x
    result%value=f
    result%iterations=min(it,nmax)
  end subroutine bfgs_minimize

  pure function outer(a,b) result(c)
    real(dp), intent(in) :: a(:),b(:)
    real(dp) :: c(size(a),size(b))
    integer :: i,j
    do j=1,size(b)
    do i=1,size(a)
    c(i,j)=a(i)*b(j)
    end do
    end do
  end function outer

  subroutine invert_matrix(a,ainv,ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    logical, intent(out), optional :: ok
    real(dp), allocatable :: aug(:,:)
    real(dp) :: piv,tmp
    integer :: n,i,j,k,imax
    logical :: good
    n=size(a,1)
    good=size(a,2)==n .and. size(ainv,1)==n .and. size(ainv,2)==n
    if(.not.good) then
    if(present(ok)) ok=.false.
    return
    end if
    allocate(aug(n,2*n))
    aug(:,1:n)=a
    aug(:,n+1:)=0.0_dp
    do i=1,n
    aug(i,n+i)=1.0_dp
    end do
    do i=1,n
      imax=i
      do k=i+1,n
      if(abs(aug(k,i))>abs(aug(imax,i))) imax=k
      end do
      if(abs(aug(imax,i))<=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))) then
      good=.false.
      exit
      end if
      if(imax/=i) then
        do j=1,2*n
        tmp=aug(i,j)
        aug(i,j)=aug(imax,j)
        aug(imax,j)=tmp
        end do
      end if
      piv=aug(i,i)
      aug(i,:)=aug(i,:)/piv
      do k=1,n
        if(k/=i) aug(k,:)=aug(k,:)-aug(k,i)*aug(i,:)
      end do
    end do
    if(good) ainv=aug(:,n+1:)
    if(present(ok)) ok=good
  end subroutine invert_matrix
end module gb2_optimizer
