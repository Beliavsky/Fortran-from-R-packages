! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_fracpoly
  use flexsurv_kinds, only : dp
  implicit none
  private
  public :: bfp, dbfp
contains
  subroutine bfp(x,powers,b)
    real(dp),intent(in)::x(:),powers(:)
    real(dp),intent(out)::b(size(x),size(powers))
    real(dp),allocatable::x1(:),x2(:)
    integer::i
    allocate(x1(size(x)),x2(size(x)))
    if(abs(powers(1))<=epsilon(1.0_dp))then;x1=log(x);else;x1=x**powers(1);end if
    b(:,1)=x1
    do i=2,size(powers)
      if(abs(powers(i)-powers(i-1))<=epsilon(1.0_dp))then
        x2=log(x)*x1
      else if(abs(powers(i))<=epsilon(1.0_dp))then
        x2=log(x)
      else
        x2=x**powers(i)
      end if
      b(:,i)=x2;x1=x2
    end do
  end subroutine bfp

  subroutine dbfp(x,powers,b)
    real(dp),intent(in)::x(:),powers(:)
    real(dp),intent(out)::b(size(x),size(powers))
    real(dp),allocatable::x1(:),x2(:),dx1(:),dx2(:)
    integer::i
    allocate(x1(size(x)),x2(size(x)),dx1(size(x)),dx2(size(x)))
    if(abs(powers(1))<=epsilon(1.0_dp))then
      x1=log(x);dx1=1.0_dp/x
    else
      x1=x**powers(1);dx1=powers(1)*x**(powers(1)-1.0_dp)
    end if
    b(:,1)=dx1
    do i=2,size(powers)
      if(abs(powers(i)-powers(i-1))<=epsilon(1.0_dp))then
        x2=log(x)*x1;dx2=log(x)*dx1+x1/x
      else if(abs(powers(i))<=epsilon(1.0_dp))then
        x2=log(x);dx2=1.0_dp/x
      else
        x2=x**powers(i);dx2=powers(i)*x**(powers(i)-1.0_dp)
      end if
      b(:,i)=dx2;x1=x2;dx1=dx2
    end do
  end subroutine dbfp
end module flexsurv_fracpoly
