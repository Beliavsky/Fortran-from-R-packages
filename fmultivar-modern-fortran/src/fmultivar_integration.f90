! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fmultivar_integration
  use fmultivar_kinds, only : dp
  implicit none
  private
  public :: integration_result, integrate_1d, integrate2d_rule, adapt_integrate2d, adapt_integrate_nd

  type :: integration_result
    real(dp) :: value = 0.0_dp
    real(dp) :: error = 0.0_dp
    integer :: evaluations = 0
    logical :: converged = .false.
  end type integration_result

  abstract interface
    function scalar_fun1d(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_fun1d
    function scalar_fun2d(x,y) result(z)
      import dp
      real(dp), intent(in) :: x,y
      real(dp) :: z
    end function scalar_fun2d
    function scalar_fun_nd(x) result(z)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: z
    end function scalar_fun_nd
  end interface
contains
  function integrate_1d(fun,a,b,tol,max_depth) result(res)
    procedure(scalar_fun1d) :: fun
    real(dp), intent(in) :: a,b,tol
    integer, intent(in), optional :: max_depth
    type(integration_result) :: res
    integer :: depth
    real(dp) :: fa,fb,fm,s
    depth=20; if (present(max_depth)) depth=max_depth
    fa=fun(a); fb=fun(b); fm=fun(0.5_dp*(a+b))
    res%evaluations=3
    s=(b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
    call adapt_simpson(fun,a,b,fa,fm,fb,s,tol,depth,res%value,res%error,res%evaluations,res%converged)
  end function integrate_1d

  recursive subroutine adapt_simpson(fun,a,b,fa,fm,fb,s,tol,depth,val,err,neval,ok)
    procedure(scalar_fun1d) :: fun
    real(dp), intent(in) :: a,b,fa,fm,fb,s,tol
    integer, intent(in) :: depth
    real(dp), intent(out) :: val,err
    integer, intent(inout) :: neval
    logical, intent(out) :: ok
    real(dp) :: m,lm,rm,flm,frm,sl,sr,delta,v1,v2,e1,e2
    logical :: ok1,ok2
    m=0.5_dp*(a+b); lm=0.5_dp*(a+m); rm=0.5_dp*(m+b)
    flm=fun(lm); frm=fun(rm); neval=neval+2
    sl=(m-a)*(fa+4.0_dp*flm+fm)/6.0_dp
    sr=(b-m)*(fm+4.0_dp*frm+fb)/6.0_dp
    delta=sl+sr-s
    if (depth<=0 .or. abs(delta)<=15.0_dp*tol) then
      val=sl+sr+delta/15.0_dp; err=abs(delta)/15.0_dp; ok=depth>0
    else
      call adapt_simpson(fun,a,m,fa,flm,fm,sl,0.5_dp*tol,depth-1,v1,e1,neval,ok1)
      call adapt_simpson(fun,m,b,fm,frm,fb,sr,0.5_dp*tol,depth-1,v2,e2,neval,ok2)
      val=v1+v2; err=e1+e2; ok=ok1.and.ok2
    end if
  end subroutine adapt_simpson

  function integrate2d_rule(fun,error_target) result(res)
    procedure(scalar_fun2d) :: fun
    real(dp), intent(in), optional :: error_target
    type(integration_result) :: res
    real(dp) :: target,h,x0,y0
    integer :: n,blocks,i,j
    target=1.0e-5_dp; if(present(error_target)) target=error_target
    h=sqrt(sqrt(max(target,epsilon(1.0_dp))))
    n=ceiling(1.0_dp/h+1.0_dp)
    blocks=ceiling(log(real(n+1,dp))/log(2.0_dp))
    n=2**blocks-1
    if(mod(n,2)==0)n=n+1
    h=1.0_dp/real(n-1,dp)
    res%value=0.0_dp
    do i=2,n-1,2
      x0=real(i-1,dp)*h
      do j=2,n-1,2
        y0=real(j-1,dp)*h
        res%value=res%value+local_nine(fun,x0,y0,h)
        res%evaluations=res%evaluations+9
      end do
    end do
    res%value=4.0_dp*h*h*res%value
    res%error=h**4; res%converged=.true.
  end function integrate2d_rule

  function local_nine(fun,x,y,h) result(v)
    procedure(scalar_fun2d) :: fun
    real(dp),intent(in)::x,y,h
    real(dp)::v
    v=(16.0_dp*fun(x,y)+fun(x-h,y-h)+fun(x-h,y+h)+fun(x+h,y-h)+fun(x+h,y+h)+ &
       4.0_dp*(fun(x-h,y)+fun(x+h,y)+fun(x,y-h)+fun(x,y+h)))/36.0_dp
  end function local_nine

  function adapt_integrate2d(fun,lower,upper,tol,max_depth) result(res)
    procedure(scalar_fun2d) :: fun
    real(dp), intent(in) :: lower(2),upper(2),tol
    integer,intent(in),optional::max_depth
    type(integration_result)::res
    integer::depth
    depth=12;if(present(max_depth))depth=max_depth
    call adapt_rect(fun,lower(1),upper(1),lower(2),upper(2),tol,depth,res%value,res%error,res%evaluations,res%converged)
  end function adapt_integrate2d

  recursive subroutine adapt_rect(fun,ax,bx,ay,by,tol,depth,val,err,neval,ok)
    procedure(scalar_fun2d)::fun
    real(dp),intent(in)::ax,bx,ay,by,tol
    integer,intent(in)::depth
    real(dp),intent(out)::val,err
    integer,intent(inout)::neval
    logical,intent(out)::ok
    real(dp)::coarse,fine,mx,my,v(4),e(4)
    logical::oks(4)
    coarse=gauss5_rect(fun,ax,bx,ay,by); neval=neval+25
    mx=0.5_dp*(ax+bx);my=0.5_dp*(ay+by)
    fine=gauss5_rect(fun,ax,mx,ay,my)+gauss5_rect(fun,mx,bx,ay,my)+ &
         gauss5_rect(fun,ax,mx,my,by)+gauss5_rect(fun,mx,bx,my,by)
    neval=neval+100
    if(depth<=0 .or. abs(fine-coarse)<=tol)then
      val=fine;err=abs(fine-coarse);ok=depth>0
    else
      call adapt_rect(fun,ax,mx,ay,my,tol/4.0_dp,depth-1,v(1),e(1),neval,oks(1))
      call adapt_rect(fun,mx,bx,ay,my,tol/4.0_dp,depth-1,v(2),e(2),neval,oks(2))
      call adapt_rect(fun,ax,mx,my,by,tol/4.0_dp,depth-1,v(3),e(3),neval,oks(3))
      call adapt_rect(fun,mx,bx,my,by,tol/4.0_dp,depth-1,v(4),e(4),neval,oks(4))
      val=sum(v);err=sum(e);ok=all(oks)
    end if
  end subroutine adapt_rect

  function gauss5_rect(fun,ax,bx,ay,by) result(v)
    procedure(scalar_fun2d)::fun
    real(dp),intent(in)::ax,bx,ay,by
    real(dp)::v,x,y,hx,hy,cx,cy
    real(dp),parameter::node(5)=[-0.906179845938664_dp,-0.538469310105683_dp,0.0_dp,0.538469310105683_dp,0.906179845938664_dp]
    real(dp), parameter :: weight(5) = [ &
      0.236926885056189_dp, 0.478628670499366_dp, 0.568888888888889_dp, &
      0.478628670499366_dp, 0.236926885056189_dp ]
    integer::i,j
    hx=0.5_dp*(bx-ax);cx=0.5_dp*(bx+ax);hy=0.5_dp*(by-ay);cy=0.5_dp*(by+ay)
    v=0.0_dp
    do i = 1, 5
      x = cx + hx*node(i)
      do j = 1, 5
        y = cy + hy*node(j)
        v = v + weight(i)*weight(j)*fun(x,y)
      end do
    end do
    v=v*hx*hy
  end function gauss5_rect

  function adapt_integrate_nd(fun,lower,upper,tol,max_points) result(res)
    procedure(scalar_fun_nd)::fun
    real(dp),intent(in)::lower(:),upper(:),tol
    integer,intent(in),optional::max_points
    type(integration_result)::res
    integer::nmax,n1,n2,i,d
    real(dp)::s1,s2,vol
    real(dp),allocatable::x(:)
    nmax=262144;if(present(max_points))nmax=max_points
    n1=min(4096,max(256,nmax/4));n2=min(nmax,4*n1)
    d=size(lower)
    if(d<1 .or. d>20)error stop 'adapt_integrate_nd supports dimensions 1 through 20'
    if(size(upper)/=d)error stop 'Integration bounds have inconsistent dimensions'
    allocate(x(d));vol=product(upper-lower)
    s1=0.0_dp
    do i=1,n1
      call halton_point(i,d,x);x=lower+(upper-lower)*x;s1=s1+fun(x)
    end do
    s1=vol*s1/real(n1,dp)
    s2=0.0_dp
    do i=1,n2
      call halton_point(i,d,x);x=lower+(upper-lower)*x;s2=s2+fun(x)
    end do
    s2=vol*s2/real(n2,dp)
    res%value=s2;res%error=abs(s2-s1);res%evaluations=n1+n2;res%converged=res%error<=tol
  end function adapt_integrate_nd

  subroutine halton_point(index,d,x)
    integer,intent(in)::index,d
    real(dp),intent(out)::x(d)
    integer,parameter::primes(20)=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71]
    integer::j,n,base,digit
    real(dp)::f,r
    do j=1,d
      base=primes(j);n=index;f=1.0_dp;r=0.0_dp
      do while(n>0);f=f/real(base,dp);digit=mod(n,base);r=r+f*real(digit,dp);n=n/base;end do
      x(j)=r
    end do
  end subroutine halton_point
end module fmultivar_integration
