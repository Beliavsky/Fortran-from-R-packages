! SPDX-License-Identifier: BSD-2-Clause
module smoof_cec09
  use smoof_kinds, only : dp, pi
  implicit none
  private
  public :: uf1, uf2, uf3, uf4, uf5, uf6, uf7, uf8, uf9, uf10
contains

  pure subroutine uf1(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(2)
    integer::j,c1,c2,n;real(dp)::s1,s2,y
    n=size(x);s1=0;s2=0;c1=0;c2=0
    do j=2,n
      y=x(j)-sin(6.0_dp*pi*x(1)+real(j,dp)*pi/real(n,dp));y=y*y
      if(mod(j,2)==0)then;s2=s2+y;c2=c2+1;else;s1=s1+y;c1=c1+1;end if
    end do
    f(1)=x(1)+2.0_dp*s1/real(c1,dp)
    f(2)=1.0_dp-sqrt(x(1))+2.0_dp*s2/real(c2,dp)
  end subroutine uf1

  pure subroutine uf2(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(2)
    integer::j,c1,c2,n;real(dp)::s1,s2,y,t
    n=size(x);s1=0;s2=0;c1=0;c2=0
    do j=2,n
      t=6.0_dp*pi*x(1)+real(j,dp)*pi/real(n,dp)
      if(mod(j,2)==0)then
        y=x(j)-0.3_dp*x(1)*(x(1)*cos(4.0_dp*t)+2.0_dp)*cos(t);s2=s2+y*y;c2=c2+1
      else
        y=x(j)-0.3_dp*x(1)*(x(1)*cos(4.0_dp*t)+2.0_dp)*sin(t);s1=s1+y*y;c1=c1+1
      end if
    end do
    f(1)=x(1)+2.0_dp*s1/real(c1,dp)
    f(2)=1.0_dp-sqrt(x(1))+2.0_dp*s2/real(c2,dp)
  end subroutine uf2

  pure subroutine uf3(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(2)
    integer::j,c1,c2,n;real(dp)::s1,s2,p1,p2,y,p
    n=size(x);s1=0;s2=0;p1=1;p2=1;c1=0;c2=0
    do j=2,n
      y=x(j)-x(1)**(0.5_dp*(1.0_dp+3.0_dp*real(j-2,dp)/real(n-2,dp)))
      p=cos(20.0_dp*y*pi/sqrt(real(j,dp)))
      if(mod(j,2)==0)then;s2=s2+y*y;p2=p2*p;c2=c2+1;else;s1=s1+y*y;p1=p1*p;c1=c1+1;end if
    end do
    f(1)=x(1)+2.0_dp*(4.0_dp*s1-2.0_dp*p1+2.0_dp)/real(c1,dp)
    f(2)=1.0_dp-sqrt(x(1))+2.0_dp*(4.0_dp*s2-2.0_dp*p2+2.0_dp)/real(c2,dp)
  end subroutine uf3

  pure subroutine uf4(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(2)
    integer::j,c1,c2,n;real(dp)::s1,s2,y,h
    n=size(x);s1=0;s2=0;c1=0;c2=0
    do j=2,n
      y=x(j)-sin(6.0_dp*pi*x(1)+real(j,dp)*pi/real(n,dp))
      h=abs(y)/(1.0_dp+exp(2.0_dp*abs(y)))
      if(mod(j,2)==0)then;s2=s2+h;c2=c2+1;else;s1=s1+h;c1=c1+1;end if
    end do
    f(1)=x(1)+2.0_dp*s1/real(c1,dp)
    f(2)=1.0_dp-x(1)**2+2.0_dp*s2/real(c2,dp)
  end subroutine uf4

  pure subroutine uf5(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(2)
    integer::j,c1,c2,n;real(dp)::s1,s2,y,h
    n=size(x);s1=0;s2=0;c1=0;c2=0
    do j=2,n
      y=x(j)-sin(6.0_dp*pi*x(1));h=2.0_dp*y*y-cos(4.0_dp*pi*y)+1.0_dp
      if(mod(j,2)==0)then;s2=s2+h;c2=c2+1;else;s1=s1+h;c1=c1+1;end if
    end do
    h=0.5_dp*(0.05_dp+0.1_dp)*abs(sin(20.0_dp*pi*x(1)))
    f(1)=x(1)+h+2.0_dp*s1/real(c1,dp)
    f(2)=1.0_dp-x(1)+h+2.0_dp*s2/real(c2,dp)
  end subroutine uf5

  pure subroutine uf6(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(2)
    integer::j,c1,c2,n;real(dp)::s1,s2,y,h
    n=size(x);s1=0;s2=0;c1=0;c2=0
    do j=2,n
      y=x(j)-sin(6.0_dp*pi*x(1)+real(j,dp)*pi/real(n,dp))
      if(mod(j,2)==0)then;s2=s2+y*y;c2=c2+1;else;s1=s1+y*y;c1=c1+1;end if
    end do
    h=0.5_dp*(0.25_dp+0.1_dp)*sin(4.0_dp*pi*x(1));h=max(h,0.0_dp)
    f(1)=x(1)+h+2.0_dp*s1/real(c1,dp)
    f(2)=1.0_dp-x(1)+h+2.0_dp*s2/real(c2,dp)
  end subroutine uf6

  pure subroutine uf7(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(2)
    integer::j,c1,c2,n;real(dp)::s1,s2,y,z
    n=size(x);s1=0;s2=0;c1=0;c2=0
    do j=2,n
      y=x(j)-sin(6.0_dp*pi*x(1)+real(j,dp)*pi/real(n,dp))
      if(mod(j,2)==0)then;s2=s2+y*y;c2=c2+1;else;s1=s1+y*y;c1=c1+1;end if
    end do
    z=x(1)**0.2_dp
    f(1)=z+2.0_dp*s1/real(c1,dp);f(2)=1.0_dp-z+2.0_dp*s2/real(c2,dp)
  end subroutine uf7

  pure subroutine uf8(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(3)
    integer::j,c1,c2,c3,n;real(dp)::s1,s2,s3,y
    n=size(x);s1=0;s2=0;s3=0;c1=0;c2=0;c3=0
    do j=3,n
      y=x(j)-2.0_dp*x(2)*sin(2.0_dp*pi*x(1)+real(j,dp)*pi/real(n,dp))
      select case(mod(j,3));case(1);s1=s1+y*y;c1=c1+1;case(2);s2=s2+y*y;c2=c2+1;case default;s3=s3+y*y;c3=c3+1;end select
    end do
    f(1)=cos(0.5_dp*pi*x(1))*cos(0.5_dp*pi*x(2))+2.0_dp*s1/real(c1,dp)
    f(2)=cos(0.5_dp*pi*x(1))*sin(0.5_dp*pi*x(2))+2.0_dp*s2/real(c2,dp)
    f(3)=sin(0.5_dp*pi*x(1))+2.0_dp*s3/real(c3,dp)
  end subroutine uf8

  pure subroutine uf9(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(3)
    integer::j,c1,c2,c3,n;real(dp)::s1,s2,s3,y
    n=size(x);s1=0;s2=0;s3=0;c1=0;c2=0;c3=0
    do j=3,n
      y=x(j)-2.0_dp*x(2)*sin(2.0_dp*pi*x(1)+real(j,dp)*pi/real(n,dp))
      select case(mod(j,3));case(1);s1=s1+y*y;c1=c1+1;case(2);s2=s2+y*y;c2=c2+1;case default;s3=s3+y*y;c3=c3+1;end select
    end do
    y=max((0.5_dp+0.1_dp)*(1.0_dp-4.0_dp*(2.0_dp*x(1)-1.0_dp)**2),0.0_dp)
    f(1)=0.5_dp*(y+2.0_dp*x(1))*x(2)+2.0_dp*s1/real(c1,dp)
    f(2)=0.5_dp*(y-2.0_dp*x(1)+2.0_dp)*x(2)+2.0_dp*s2/real(c2,dp)
    f(3)=1.0_dp-x(2)+2.0_dp*s3/real(c3,dp)
  end subroutine uf9

  pure subroutine uf10(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(3)
    integer::j,c1,c2,c3,n;real(dp)::s1,s2,s3,y,h
    n=size(x);s1=0;s2=0;s3=0;c1=0;c2=0;c3=0
    do j=3,n
      y=x(j)-2.0_dp*x(2)*sin(2.0_dp*pi*x(1)+real(j,dp)*pi/real(n,dp))
      h=4.0_dp*y*y-cos(8.0_dp*pi*y)+1.0_dp
      select case(mod(j,3));case(1);s1=s1+h;c1=c1+1;case(2);s2=s2+h;c2=c2+1;case default;s3=s3+h;c3=c3+1;end select
    end do
    f(1)=cos(0.5_dp*pi*x(1))*cos(0.5_dp*pi*x(2))+2.0_dp*s1/real(c1,dp)
    f(2)=cos(0.5_dp*pi*x(1))*sin(0.5_dp*pi*x(2))+2.0_dp*s2/real(c2,dp)
    f(3)=sin(0.5_dp*pi*x(1))+2.0_dp*s3/real(c3,dp)
  end subroutine uf10
end module smoof_cec09
