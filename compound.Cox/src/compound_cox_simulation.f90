! SPDX-License-Identifier: GPL-2.0-only
module compound_cox_simulation
  use compound_cox_kinds, only : dp
  implicit none
  private
  public :: x_pathway, x_tag
contains
  subroutine x_pathway(n,p,q1,q2,x,rho1,rho2)
    integer,intent(in)::n,p,q1,q2
    real(dp),allocatable,intent(out)::x(:,:)
    real(dp),intent(in),optional::rho1,rho2
    real(dp)::rh1,rh2,r1,r2,sd0,sd1,sd2,a,u
    integer::i,j
    if(q1+q2>p)error stop 'x_pathway: q1+q2 > p'
    rh1=0.5_dp
    if(present(rho1))rh1=rho1
    rh2=0.5_dp
    if(present(rho2))rh2=rho2
    r1=1.0_dp/(1.0_dp+sqrt((1.0_dp-rh1)/rh1))
    r2=1.0_dp/(1.0_dp+sqrt((1.0_dp-rh2)/rh2))
    sd0=sqrt(3.0_dp/4.0_dp)
    sd1=sqrt(0.75_dp*(2*r1*r1-2*r1+1))
    sd2=sqrt(0.75_dp*(2*r2*r2-2*r2+1))
    allocate(x(n,p))
    x=0
    do i=1,n
      if(q1>0)then
      call random_number(u)
      a=-1.5_dp*r1+3.0_dp*r1*u
      do j=1,q1
      call random_number(u)
      x(i,j)=(-1.5_dp*(1-r1)+3.0_dp*(1-r1)*u+a)/sd1
      end do
      end if
      if(q2>0)then
      call random_number(u)
      a=-1.5_dp*r2+3.0_dp*r2*u
      do j=q1+1,q1+q2
      call random_number(u)
      x(i,j)=(-1.5_dp*(1-r2)+3.0_dp*(1-r2)*u+a)/sd2
      end do
      end if
      do j=q1+q2+1,p
      call random_number(u)
      x(i,j)=(-1.5_dp+3.0_dp*u)/sd0
      end do
    end do
  end subroutine x_pathway

  subroutine x_tag(n,p,q,x,s)
    integer,intent(in)::n,p,q
    real(dp),allocatable,intent(out)::x(:,:)
    integer,intent(in),optional::s
    integer::ss,i,j,k,col
    real(dp)::a,u
    ss=1
    if(present(s))ss=s
    if(q+q*ss>p)error stop 'x_tag: q+q*s > p'
    allocate(x(n,p))
    x=0
    do i=1,n
      do j=1,q
        call random_number(u)
        a=-0.75_dp+1.5_dp*u
        col=j
        call random_number(u)
        x(i,col)=(-0.75_dp+1.5_dp*u+a)/0.612_dp
        do k=1,ss
        col=q+ss*(j-1)+k
        call random_number(u)
        x(i,col)=(-0.75_dp+1.5_dp*u+a)/0.612_dp
        end do
      end do
      do col=q+q*ss+1,p
      call random_number(u)
      x(i,col)=(-1.5_dp+3.0_dp*u)/0.866_dp
      end do
    end do
  end subroutine x_tag
end module compound_cox_simulation
