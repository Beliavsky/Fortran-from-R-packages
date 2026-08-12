! SPDX-License-Identifier: GPL-3.0-or-later
! Modern Fortran translation of computational code from R package adagio 0.9.2.
module adagio_maxquad
  use adagio_kinds, only : dp
  implicit none
  private
  public :: maxquad_problem, make_maxquad

  type :: maxquad_problem
     integer :: n=0,m=0
     real(dp),allocatable :: a(:,:,:)
     real(dp),allocatable :: b(:,:)
   contains
     procedure :: value => maxquad_value
     procedure :: gradient => maxquad_gradient
  end type

contains

  function make_maxquad(n,m) result(prob)
    integer,intent(in)::n,m
    type(maxquad_problem)::prob
    integer::l,i,j
    prob%n=n;prob%m=m;allocate(prob%a(n,n,m),prob%b(n,m));prob%a=0
    do l=1,m
       do i=1,n
          do j=i+1,n
             prob%a(i,j,l)=exp(real(i,dp)/real(j,dp))*cos(real(i*j,dp))*sin(real(l,dp))
             prob%a(j,i,l)=prob%a(i,j,l)
          end do
          prob%b(i,l)=exp(real(i,dp)/real(l,dp))*sin(real(i*l,dp))
       end do
       do i=1,n
          prob%a(i,i,l)=real(i,dp)*abs(sin(real(l,dp)))/real(n,dp)+sum(abs(prob%a(i,:,l)))
       end do
    end do
  end function make_maxquad

  function maxquad_value(self,x) result(f)
    class(maxquad_problem),intent(in)::self
    real(dp),intent(in)::x(:)
    real(dp)::f,d
    integer::l
    f=huge(1.0_dp)
    do l=1,self%m
       d=dot_product(x,matmul(self%a(:,:,l),x))-dot_product(self%b(:,l),x)
       if(l==1 .or. d>f) f=d
    end do
  end function maxquad_value

  function maxquad_gradient(self,x) result(g)
    class(maxquad_problem),intent(in)::self
    real(dp),intent(in)::x(:)
    real(dp)::g(size(x)),f,d
    integer::l,k
    f=-huge(1.0_dp);k=1
    do l=1,self%m
       d=dot_product(x,matmul(self%a(:,:,l),x))-dot_product(self%b(:,l),x)
       if(d>f)then;f=d;k=l;end if
    end do
    g=2.0_dp*matmul(self%a(:,:,k),x)-self%b(:,k)
  end function maxquad_gradient
end module adagio_maxquad
