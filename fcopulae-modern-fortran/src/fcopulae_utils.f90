! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fcopulae_utils
  use fcopulae_kinds, only : dp, pi
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: clamp01, pseudo_observations, kendall_tau_sample, spearman_rho_sample
  public :: numerical_first, numerical_second, numerical_mixed_density
  public :: debye_function, safe_log, is_finite_real
  abstract interface
    function scalar_function(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_function
    function bivariate_function(x,y) result(z)
      import dp
      real(dp), intent(in) :: x,y
      real(dp) :: z
    end function bivariate_function
  end interface
contains
  elemental function clamp01(x,eps) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::eps
    real(dp)::y,e
    e=1.0e-12_dp;if(present(eps))e=eps
    y=min(1.0_dp-e,max(e,x))
  end function clamp01

  elemental function safe_log(x) result(y)
    real(dp),intent(in)::x
    real(dp)::y
    y=log(max(x,tiny(1.0_dp)))
  end function safe_log

  elemental logical function is_finite_real(x)
    real(dp),intent(in)::x
    is_finite_real=ieee_is_finite(x)
  end function is_finite_real

  subroutine pseudo_observations(x,u)
    real(dp),intent(in)::x(:)
    real(dp),allocatable,intent(out)::u(:)
    integer::n,i,j,less,equal
    n=size(x);allocate(u(n))
    do i=1,n
      less=0;equal=0
      do j=1,n
        if(x(j)<x(i))less=less+1
        if(abs(x(j)-x(i))<=epsilon(1.0_dp)*max(1.0_dp,abs(x(i)),abs(x(j))))equal=equal+1
      end do
      u(i)=(real(less,dp)+0.5_dp*real(equal+1,dp))/real(n+1,dp)
    end do
  end subroutine pseudo_observations

  function kendall_tau_sample(x,y) result(tau)
    real(dp),intent(in)::x(:),y(:)
    real(dp)::tau
    integer::i,j,n
    real(dp)::s,dx,dy
    n=min(size(x),size(y));s=0.0_dp
    do i=1,n-1
      do j=i+1,n
        dx=x(i)-x(j);dy=y(i)-y(j)
        if(dx*dy>0.0_dp)s=s+1.0_dp
        if(dx*dy<0.0_dp)s=s-1.0_dp
      end do
    end do
    if(n<2)then;tau=0.0_dp;else;tau=2.0_dp*s/real(n*(n-1),dp);end if
  end function kendall_tau_sample

  function spearman_rho_sample(x,y) result(rho)
    real(dp),intent(in)::x(:),y(:)
    real(dp)::rho
    real(dp),allocatable::rx(:),ry(:)
    real(dp)::mx,my,sx,sy
    integer::n
    call pseudo_observations(x,rx);call pseudo_observations(y,ry)
    n=min(size(rx),size(ry));mx=sum(rx(1:n))/real(n,dp);my=sum(ry(1:n))/real(n,dp)
    sx=sqrt(sum((rx(1:n)-mx)**2));sy=sqrt(sum((ry(1:n)-my)**2))
    if(sx<=0.0_dp .or. sy<=0.0_dp)then;rho=0.0_dp
    else;rho=sum((rx(1:n)-mx)*(ry(1:n)-my))/(sx*sy);end if
  end function spearman_rho_sample

  function numerical_first(fun,x,lower,upper) result(d)
    procedure(scalar_function)::fun
    real(dp),intent(in)::x,lower,upper
    real(dp)::d,h,xm,xp
    h=max(1.0e-6_dp,1.0e-5_dp*max(1.0_dp,abs(x)))
    xm=max(lower,x-h);xp=min(upper,x+h)
    if(xp<=xm)then;d=0.0_dp;else;d=(fun(xp)-fun(xm))/(xp-xm);end if
  end function numerical_first

  function numerical_second(fun,x,lower,upper) result(d2)
    procedure(scalar_function)::fun
    real(dp),intent(in)::x,lower,upper
    real(dp)::d2,h,xm,xp,h1,h2,f0,fm,fp
    h=max(2.0e-5_dp,5.0e-5_dp*max(1.0_dp,abs(x)))
    xm=max(lower,x-h);xp=min(upper,x+h);h1=x-xm;h2=xp-x
    if(h1<=0.0_dp .or. h2<=0.0_dp)then;d2=0.0_dp;return;end if
    fm=fun(xm);f0=fun(x);fp=fun(xp)
    d2=2.0_dp*(fm/(h1*(h1+h2))-f0/(h1*h2)+fp/(h2*(h1+h2)))
  end function numerical_second

  function numerical_mixed_density(fun,u,v) result(d)
    procedure(bivariate_function)::fun
    real(dp),intent(in)::u,v
    real(dp)::d,h,um,up,vm,vp
    h=2.0e-5_dp
    um=max(1.0e-10_dp,u-h);up=min(1.0_dp-1.0e-10_dp,u+h)
    vm=max(1.0e-10_dp,v-h);vp=min(1.0_dp-1.0e-10_dp,v+h)
    d=(fun(up,vp)-fun(up,vm)-fun(um,vp)+fun(um,vm))/((up-um)*(vp-vm))
    d=max(0.0_dp,d)
  end function numerical_mixed_density

  function debye_function(x,k) result(d)
    real(dp),intent(in)::x
    integer,intent(in),optional::k
    real(dp)::d,u,h,t,s
    integer::kk,n,i
    kk=1;if(present(k))kk=k
    if(kk<=0)then;d=0.0_dp;return;end if
    if(abs(x)<1.0e-12_dp)then;d=1.0_dp;return;end if
    u=abs(x);n=2000;if(mod(n,2)/=0)n=n+1;h=u/real(n,dp);s=0.0_dp
    do i=0,n
      t=real(i,dp)*h
      if(i==0)then
        if(kk==1)then;s=s+1.0_dp;else;s=s+0.0_dp;end if
      else if(i==n)then
        s=s+t**kk/(exp(t)-1.0_dp)
      else if(mod(i,2)==0)then
        s=s+2.0_dp*t**kk/(exp(t)-1.0_dp)
      else
        s=s+4.0_dp*t**kk/(exp(t)-1.0_dp)
      end if
    end do
    d=real(kk,dp)*(h*s/3.0_dp)/u**kk
    if(x<0.0_dp)d=d+real(kk,dp)*u/real(kk+1,dp)
  end function debye_function
end module fcopulae_utils
