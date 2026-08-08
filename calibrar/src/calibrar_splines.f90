! SPDX-License-Identifier: GPL-2.0-only
module calibrar_splines
  use calibrar_kinds, only : dp
  implicit none
  private
  public :: spline_result, spline_par, cubic_spline_eval

  type :: spline_result
    real(dp), allocatable :: time(:)
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: knots(:)
    real(dp), allocatable :: par(:)
    logical :: periodic = .false.
  end type spline_result

contains

  subroutine spline_par(par, n, result, knots, periodic, period)
    real(dp), intent(in) :: par(:)
    integer, intent(in) :: n
    type(spline_result), intent(out) :: result
    real(dp), intent(in), optional :: knots(:)
    logical, intent(in), optional :: periodic
    real(dp), intent(in), optional :: period
    logical :: per
    real(dp), allocatable :: k0(:), p0(:), kk(:), pp(:)
    real(dp) :: den, perval
    integer :: i, m
    per=.false.; if(present(periodic)) per=periodic
    m=size(par)
    if(m<2 .or. n<1) error stop "spline_par: invalid dimensions"
    allocate(k0(m), p0(m)); p0=par
    if(present(knots)) then
      if(size(knots)/=m) error stop "spline_par: knots and par size mismatch"
      k0=knots
      do i=2,m
        if(k0(i)<=k0(i-1)) error stop "spline_par: knots must be strictly increasing"
      end do
      if(per) then
        if(.not.present(period)) error stop "spline_par: period required with explicit periodic knots"
        perval=period
        if(perval<=0.0_dp) error stop "spline_par: period must be positive"
        k0=modulo(k0/perval,1.0_dp)
        call sort_pairs(k0,p0)
      else
        den=maxval(k0)-minval(k0)
        if(den<=0.0_dp) error stop "spline_par: degenerate knots"
        k0=(k0-minval(k0))/den
      end if
    else
      do i=1,m
        k0(i)=real(i-1,dp)/real(m-1,dp)
      end do
    end if
    if(per) then
      allocate(kk(3*m),pp(3*m))
      kk(1:m)=k0-1.0_dp; kk(m+1:2*m)=k0; kk(2*m+1:3*m)=k0+1.0_dp
      pp(1:m)=p0; pp(m+1:2*m)=p0; pp(2*m+1:3*m)=p0
    else
      allocate(kk(m),pp(m)); kk=k0; pp=p0
    end if
    allocate(result%time(n),result%x(n),result%knots(m),result%par(m))
    do i=1,n
      result%time(i)=(real(i,dp)-0.5_dp)/real(n,dp)
    end do
    call cubic_spline_eval(kk,pp,result%time,result%x)
    result%knots=k0; result%par=p0; result%periodic=per
  end subroutine spline_par

  subroutine cubic_spline_eval(x, y, xq, yq)
    real(dp), intent(in) :: x(:), y(:), xq(:)
    real(dp), intent(out) :: yq(:)
    real(dp), allocatable :: a(:,:), rhs(:), m2(:)
    real(dp) :: h1,h2,t,h
    integer :: n,i,j,k
    n=size(x)
    if(size(y)/=n .or. size(yq)/=size(xq)) error stop "cubic_spline_eval: size mismatch"
    if(n==2) then
      do i=1,size(xq)
        yq(i)=y(1)+(y(2)-y(1))*(xq(i)-x(1))/(x(2)-x(1))
      end do
      return
    end if
    allocate(a(n,n),rhs(n),m2(n)); a=0.0_dp; rhs=0.0_dp
    if(n>=4) then
      h1=x(2)-x(1); h2=x(3)-x(2)
      a(1,1)=-h2; a(1,2)=h1+h2; a(1,3)=-h1
      h1=x(n-1)-x(n-2); h2=x(n)-x(n-1)
      a(n,n-2)=-h2; a(n,n-1)=h1+h2; a(n,n)=-h1
    else
      a(1,1)=1.0_dp; a(n,n)=1.0_dp
    end if
    do i=2,n-1
      h1=x(i)-x(i-1); h2=x(i+1)-x(i)
      a(i,i-1)=h1; a(i,i)=2.0_dp*(h1+h2); a(i,i+1)=h2
      rhs(i)=6.0_dp*((y(i+1)-y(i))/h2-(y(i)-y(i-1))/h1)
    end do
    call solve_dense(a,rhs,m2)
    do i=1,size(xq)
      if(xq(i)<=x(1)) then
        j=1
      else if(xq(i)>=x(n)) then
        j=n-1
      else
        j=1
        do k=1,n-1
          if(xq(i)>=x(k) .and. xq(i)<=x(k+1)) then
            j=k; exit
          end if
        end do
      end if
      h=x(j+1)-x(j); t=(xq(i)-x(j))/h
      yq(i)=(1.0_dp-t)*y(j)+t*y(j+1) + &
        (((1.0_dp-t)**3-(1.0_dp-t))*m2(j)+(t**3-t)*m2(j+1))*h*h/6.0_dp
    end do
  end subroutine cubic_spline_eval

  subroutine solve_dense(a,b,x)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    real(dp), allocatable :: aa(:,:), bb(:), rowtmp(:)
    real(dp) :: fac,tmp
    integer :: n,i,k,p
    n=size(b); allocate(aa(n,n),bb(n),rowtmp(n)); aa=a; bb=b
    do k=1,n-1
      p=k
      do i=k+1,n
        if(abs(aa(i,k))>abs(aa(p,k))) p=i
      end do
      if(abs(aa(p,k))<tiny(1.0_dp)) error stop "solve_dense: singular system"
      if(p/=k) then
        rowtmp=aa(k,:); aa(k,:)=aa(p,:); aa(p,:)=rowtmp
        tmp=bb(k); bb(k)=bb(p); bb(p)=tmp
      end if
      do i=k+1,n
        fac=aa(i,k)/aa(k,k)
        aa(i,k:n)=aa(i,k:n)-fac*aa(k,k:n)
        bb(i)=bb(i)-fac*bb(k)
      end do
    end do
    x(n)=bb(n)/aa(n,n)
    do i=n-1,1,-1
      x(i)=(bb(i)-sum(aa(i,i+1:n)*x(i+1:n)))/aa(i,i)
    end do
  end subroutine solve_dense

  subroutine sort_pairs(x,y)
    real(dp), intent(inout) :: x(:),y(:)
    real(dp) :: tx,ty
    integer :: i,j
    do i=2,size(x)
      tx=x(i); ty=y(i); j=i-1
      do while(j>=1)
        if(x(j)<=tx) exit
        x(j+1)=x(j); y(j+1)=y(j); j=j-1
      end do
      x(j+1)=tx; y(j+1)=ty
    end do
  end subroutine sort_pairs
end module calibrar_splines
