! SPDX-License-Identifier: BSD-2-Clause
module smoof_ed
  use smoof_kinds, only : dp, pi
  implicit none
  private
  public :: ed1, ed2
contains
  pure subroutine ed1(x,m,gamma,theta,f)
    real(dp),intent(in)::x(:),gamma,theta(:);integer,intent(in)::m
    real(dp),intent(out)::f(m)
    real(dp)::p(m),r
    integer::i,j
    p(1)=cos(theta(1))**(2.0_dp/gamma)
    do i=2,m-1
      p(i)=1.0_dp
      do j=1,i-1;p(i)=p(i)*sin(theta(j));end do
      p(i)=(p(i)*cos(theta(i)))**(2.0_dp/gamma)
    end do
    p(m)=1.0_dp
    do j=1,m-1;p(m)=p(m)*sin(theta(j));end do
    p(m)=p(m)**(2.0_dp/gamma)
    r=sqrt(sum(x(m:)**2));f=p/(r+1.0_dp)
  end subroutine ed1
  pure subroutine ed2(x,m,gamma,theta,f)
    real(dp),intent(in)::x(:),gamma,theta(:);integer,intent(in)::m
    real(dp),intent(out)::f(m)
    real(dp)::p(m),r,g
    integer::i,j
    p(1)=cos(theta(1))**(2.0_dp/gamma)
    do i=2,m-1
      p(i)=1.0_dp
      do j=1,i-1;p(i)=p(i)*sin(theta(j));end do
      p(i)=(p(i)*cos(theta(i)))**(2.0_dp/gamma)
    end do
    p(m)=1.0_dp
    do j=1,m-1;p(m)=p(m)*sin(theta(j));end do
    p(m)=p(m)**(2.0_dp/gamma)
    r=sqrt(sum(x(m:)**2))
    g=0.051373_dp+(r-0.051373_dp)+0.5_dp+0.5_dp*cos(2.0_dp*pi*(r-0.051373_dp)+pi)
    f=p/(g+1.0_dp)
  end subroutine ed2
end module smoof_ed
