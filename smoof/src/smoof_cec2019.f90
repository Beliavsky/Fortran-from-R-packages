! SPDX-License-Identifier: BSD-2-Clause
module smoof_cec2019
  use smoof_kinds, only : dp, pi
  implicit none
  private
  public :: sympart_simple, sympart_rotated, omni_test
  public :: mmf1, mmf1z, mmf1e, mmf2, mmf3, mmf4, mmf5, mmf6, mmf7, mmf8
  public :: mmf9, mmf10, mmf11, mmf12, mmf13, mmf14, mmf14a, mmf15, mmf15a
contains
  pure real(dp) function sgn(x) result(r)
    real(dp),intent(in)::x
    if(x>0.0_dp)then;r=1.0_dp;else if(x<0.0_dp)then;r=-1.0_dp;else;r=0.0_dp;end if
  end function sgn

  pure subroutine sympart_simple(x,a,b,c,y)
    real(dp),intent(in)::x(:),a,b,c;real(dp),intent(out)::y(2)
    real(dp)::t1h,t2h,t1,t2,p1,p2
    t1h=sgn(x(1))*real(ceiling((abs(x(1))-(a+c/2.0_dp))/(2.0_dp*a+c)),dp)
    t2h=sgn(x(2))*real(ceiling((abs(x(2))-b/2.0_dp)/b),dp)
    t1=sgn(t1h)*min(abs(t1h),1.0_dp);t2=sgn(t2h)*min(abs(t2h),1.0_dp)
    p1=x(1)-t1*(c+2.0_dp*a);p2=x(2)-t2*b
    y(1)=(p1+a)**2+p2**2;y(2)=(p1-a)**2+p2**2
  end subroutine sympart_simple

  pure subroutine sympart_rotated(x,w,a,b,c,y)
    real(dp),intent(in)::x(:),w,a,b,c;real(dp),intent(out)::y(2)
    real(dp)::r(2)
    r(1)=cos(w)*x(1)-sin(w)*x(2);r(2)=sin(w)*x(1)+cos(w)*x(2)
    call sympart_simple(r,a,b,c,y)
  end subroutine sympart_rotated

  pure subroutine omni_test(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(2)
    y(1)=sum(sin(pi*x));y(2)=sum(cos(pi*x))
  end subroutine omni_test

  pure subroutine mmf1(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(2)
    y(1)=abs(x(1)-2.0_dp)
    y(2)=1.0_dp-sqrt(y(1))+2.0_dp*(x(2)-sin(6.0_dp*pi*y(1)+pi))**2
  end subroutine mmf1

  pure subroutine mmf1z(x,k,y)
    real(dp),intent(in)::x(:),k;real(dp),intent(out)::y(2)
    y(1)=abs(x(1)-2.0_dp)
    if(x(1)<2.0_dp)then
      y(2)=1.0_dp-sqrt(y(1))+2.0_dp*(x(2)-sin(2.0_dp*k*pi*y(1)+pi))**2
    else
      y(2)=1.0_dp-sqrt(y(1))+2.0_dp*(x(2)-sin(2.0_dp*pi*y(1)+pi))**2
    end if
  end subroutine mmf1z

  pure subroutine mmf1e(x,a,y)
    real(dp),intent(in)::x(:),a;real(dp),intent(out)::y(2)
    y(1)=abs(x(1)-2.0_dp)
    if(x(1)<2.0_dp)then
      y(2)=1.0_dp-sqrt(y(1))+2.0_dp*(x(2)-sin(6.0_dp*pi*y(1)+pi))**2
    else
      y(2)=1.0_dp-sqrt(y(1))+2.0_dp*(x(2)-a**x(1)*sin(6.0_dp*pi*y(1)+pi))**2
    end if
  end subroutine mmf1e

  pure subroutine mmf2(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(2)
    real(dp)::x2,t
    y(1)=x(1);x2=x(2);if(x2>1.0_dp)x2=x2-1.0_dp
    t=x2-sqrt(x(1));y(2)=1.0_dp-sqrt(x(1))+2.0_dp*(4.0_dp*t*t-2.0_dp*cos(20.0_dp*t*pi/sqrt(2.0_dp))+2.0_dp)
  end subroutine mmf2

  pure subroutine mmf3(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(2)
    real(dp)::t
    y(1)=x(1);t=x(2)-sqrt(x(1))
    if(x(2)>=1.0_dp)then;t=t-0.5_dp;else if(x(1)<0.25_dp .and. x(2)>0.5_dp .and. x(2)<1.0_dp)then;t=t-0.5_dp;end if
    y(2)=1.0_dp-sqrt(x(1))+2.0_dp*(4.0_dp*t*t-2.0_dp*cos(20.0_dp*t*pi/sqrt(2.0_dp))+2.0_dp)
  end subroutine mmf3

  pure subroutine mmf4(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(2)
    real(dp)::x2
    x2=x(2);if(x2>=1.0_dp)x2=x2-1.0_dp
    y(1)=abs(x(1));y(2)=1.0_dp-x(1)**2+2.0_dp*(x2-sin(pi*abs(x(1))))**2
  end subroutine mmf4

  pure subroutine mmf5(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(2)
    real(dp)::x2
    y(1)=abs(x(1)-2.0_dp);x2=x(2);if(x2>1.0_dp)x2=x2-2.0_dp
    y(2)=1.0_dp-sqrt(y(1))+2.0_dp*(x2-sin(6.0_dp*pi*y(1)+pi))**2
  end subroutine mmf5

  pure subroutine mmf6(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(2)
    real(dp)::x2
    y(1)=abs(x(1)-2.0_dp);x2=x(2);if(x2>1.0_dp)x2=x2-1.0_dp
    y(2)=1.0_dp-sqrt(y(1))+2.0_dp*(x2-sin(6.0_dp*pi*y(1)+pi))**2
  end subroutine mmf6

  pure subroutine mmf7(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(2)
    real(dp)::f1,co,inner
    f1=abs(x(1)-2.0_dp);co=24.0_dp*pi*f1+4.0_dp*pi
    inner=x(2)-(0.3_dp*f1*f1*cos(co)+0.6_dp*f1)*sin(6.0_dp*pi*f1+pi)
    y(1)=f1;y(2)=1.0_dp-sqrt(f1)+inner*inner
  end subroutine mmf7

  pure subroutine mmf8(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(2)
    real(dp)::a,b,x2
    a=abs(x(1));b=sin(a);x2=x(2);if(x2>4.0_dp)x2=x2-4.0_dp
    y(1)=b;y(2)=sqrt(1.0_dp-b*b)+2.0_dp*(x2-b-a)**2
  end subroutine mmf8

  pure subroutine mmf9(x,np,y)
    real(dp),intent(in)::x(:);integer,intent(in)::np;real(dp),intent(out)::y(2)
    real(dp)::g
    y(1)=x(1);g=2.0_dp-sin(real(np,dp)*pi*x(2))**6;y(2)=g/x(1)
  end subroutine mmf9

  pure subroutine mmf10(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(2)
    real(dp)::e1,e2,g
    y(1)=x(1);e1=exp(-((x(2)-0.2_dp)/0.004_dp)**2);e2=exp(-((x(2)-0.6_dp)/0.4_dp)**2)
    g=2.0_dp-e1-0.8_dp*e2;y(2)=g/x(1)
  end subroutine mmf10

  pure subroutine mmf11(x,np,y)
    real(dp),intent(in)::x(:);integer,intent(in)::np;real(dp),intent(out)::y(2)
    real(dp)::t1,t2,g
    y(1)=x(1);t1=exp(-2.0_dp*log10(2.0_dp)*((x(2)-0.1_dp)/0.8_dp)**2)
    t2=sin(real(np,dp)*pi*x(2))**6;g=2.0_dp-t1*t2;y(2)=g/x(1)
  end subroutine mmf11

  pure subroutine mmf12(x,np,q,y)
    real(dp),intent(in)::x(:);integer,intent(in)::np,q;real(dp),intent(out)::y(2)
    real(dp)::t1,t2,g,h
    y(1)=x(1);t1=exp(-2.0_dp*log10(2.0_dp)*((x(2)-0.1_dp)/0.8_dp)**2)
    t2=sin(real(np,dp)*pi*x(2))**6;g=2.0_dp-t1*t2
    h=1.0_dp-(x(1)/g)**2-(x(1)/g)*sin(2.0_dp*pi*real(q,dp)*x(1));y(2)=g*h
  end subroutine mmf12

  pure subroutine mmf13(x,np,y)
    real(dp),intent(in)::x(:);integer,intent(in)::np;real(dp),intent(out)::y(2)
    real(dp)::t,g
    y(1)=x(1);t=x(2)+sqrt(x(3))
    g=2.0_dp-exp(-2.0_dp*log10(2.0_dp)*((t-0.1_dp)/0.8_dp)**2)*sin(real(np,dp)*pi*t)**6
    y(2)=g/x(1)
  end subroutine mmf13

  pure subroutine mmf_shape(x,m,g,y)
    real(dp),intent(in)::x(:),g;integer,intent(in)::m;real(dp),intent(out)::y(m)
    integer::i;real(dp)::gg
    gg=g
    do i=1,m-1
      y(m-i)=gg*sin(0.5_dp*pi*x(i));gg=gg*cos(0.5_dp*pi*x(i))
    end do
    y(1)=gg
  end subroutine mmf_shape

  pure subroutine mmf14(x,m,np,y)
    real(dp),intent(in)::x(:);integer,intent(in)::m,np;real(dp),intent(out)::y(m)
    real(dp)::g
    g=3.0_dp-sin(real(np,dp)*pi*x(size(x)))**2;call mmf_shape(x,m,g,y)
  end subroutine mmf14

  pure subroutine mmf14a(x,m,np,y)
    real(dp),intent(in)::x(:);integer,intent(in)::m,np;real(dp),intent(out)::y(m)
    real(dp)::t,g
    t=x(size(x))-0.5_dp*sin(pi*x(size(x)-1))+0.5_dp/real(np,dp)
    g=3.0_dp-sin(real(np,dp)*pi*t)**2;call mmf_shape(x,m,g,y)
  end subroutine mmf14a

  pure subroutine mmf15(x,m,np,y)
    real(dp),intent(in)::x(:);integer,intent(in)::m,np;real(dp),intent(out)::y(m)
    real(dp)::t,g
    t=x(size(x));g=3.0_dp-exp(-2.0_dp*log10(2.0_dp)*((t-0.1_dp)/0.8_dp)**2)*sin(real(np,dp)*pi*t)**2
    call mmf_shape(x,m,g,y)
  end subroutine mmf15

  pure subroutine mmf15a(x,m,np,y)
    real(dp),intent(in)::x(:);integer,intent(in)::m,np;real(dp),intent(out)::y(m)
    real(dp)::t,g
    t=x(size(x))-0.5_dp*sin(pi*x(size(x)-1))+0.5_dp/real(np,dp)
    g=3.0_dp-exp(-2.0_dp*log10(2.0_dp)*((t-0.1_dp)/0.8_dp)**2)*sin(real(np,dp)*pi*t)**2
    call mmf_shape(x,m,g,y)
  end subroutine mmf15a
end module smoof_cec2019
