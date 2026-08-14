! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_benchmarks
  use gpareto_kinds, only : dp, pi
  implicit none
  private
  public :: zdt1,zdt2,zdt3,zdt4,zdt6,p1_test,p2_test,mop2,mop3,dtlz1,dtlz2,dtlz3,dtlz7,oka1
contains
  subroutine zdt1(x,y)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::g
    integer::i,n
    n=size(x,2)
    allocate(y(size(x,1),2))
    do i=1,size(x,1)
    g=1.0_dp+9.0_dp*sum(x(i,2:n))/real(n-1,dp)
    y(i,1)=x(i,1)
    y(i,2)=g*(1.0_dp-sqrt(x(i,1)/g))
    end do
  end subroutine
  subroutine zdt2(x,y)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::g
    integer::i,n
    n=size(x,2)
    allocate(y(size(x,1),2))
    do i=1,size(x,1)
    g=1.0_dp+9.0_dp*sum(x(i,2:n))/real(n-1,dp)
    y(i,1)=x(i,1)
    y(i,2)=g*(1.0_dp-(x(i,1)/g)**2)
    end do
  end subroutine
  subroutine zdt3(x,y)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::g
    integer::i,n
    n=size(x,2)
    allocate(y(size(x,1),2))
    do i=1,size(x,1)
    g=1.0_dp+9.0_dp*sum(x(i,2:n))/real(n-1,dp)
    y(i,1)=x(i,1)
    y(i,2)=g*(1.0_dp-sqrt(x(i,1)/g)-x(i,1)/g*sin(10.0_dp*pi*x(i,1)))
    end do
  end subroutine
  subroutine zdt4(x,y)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::g,t
    integer::i,j,n
    n=size(x,2)
    allocate(y(size(x,1),2))
    do i=1,size(x,1)
    g=1.0_dp+10.0_dp*real(n-1,dp)
    do j=2,n
    t=10.0_dp*x(i,j)-5.0_dp
    g=g+t*t-10.0_dp*cos(4.0_dp*pi*t)
    end do
    y(i,1)=x(i,1)
    y(i,2)=g*(1.0_dp-sqrt(x(i,1)/g))
    end do
  end subroutine
  subroutine zdt6(x,y)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::g,f1
    integer::i,n
    n=size(x,2)
    allocate(y(size(x,1),2))
    do i=1,size(x,1)
    f1=1.0_dp-exp(-4.0_dp*x(i,1))*sin(6.0_dp*pi*x(i,1))**6
    g=1.0_dp+9.0_dp*(sum(x(i,2:n))/real(n-1,dp))**0.25_dp
    y(i,:)=[f1,g*(1.0_dp-(f1/g)**2)]
    end do
  end subroutine
  subroutine p1_test(x,y)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::b1,b2
    integer::i
    allocate(y(size(x,1),2))
    do i=1,size(x,1)
    b1=15.0_dp*x(i,1)-5.0_dp
    b2=15.0_dp*x(i,2)
    y(i,1)=(b2-5.1_dp*(b1/(2*pi))**2+5.0_dp/pi*b1-6.0_dp)**2+10.0_dp*((1.0_dp-1.0_dp/(8*pi))*cos(b1)+1.0_dp)
    y(i,2)=-sqrt((10.5_dp-b1)*(b1+5.5_dp)*(b2+0.5_dp)) &
      -(b2-5.1_dp*(b1/(2*pi))**2-6.0_dp)**2/30.0_dp &
      -((1.0_dp-1.0_dp/(8*pi))*cos(b1)+1.0_dp)/3.0_dp
    end do
  end subroutine
  subroutine p2_test(x,y)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp),parameter::ap(2,2)=reshape([0.5_dp,1.5_dp,1.0_dp,2.0_dp],[2,2])
    real(dp),parameter::bp(2,2)=reshape([-2.0_dp,-1.0_dp,-1.5_dp,-0.5_dp],[2,2])
    real(dp)::a1,a2,b1,b2,u,v,f1,f2
    integer::i
    allocate(y(size(x,1),2))
    a1=ap(1,1)*sin(1.0_dp)+bp(1,1)*cos(1.0_dp)+ap(1,2)*sin(2.0_dp)+bp(1,2)*cos(2.0_dp)
    a2=ap(2,1)*sin(1.0_dp)+bp(2,1)*cos(1.0_dp)+ap(2,2)*sin(2.0_dp)+bp(2,2)*cos(2.0_dp)
    do i=1,size(x,1)
    u=2*pi*x(i,1)-pi
    v=2*pi*x(i,2)-pi
    b1=ap(1,1)*sin(u)+bp(1,1)*cos(u)+ap(1,2)*sin(v)+bp(1,2)*cos(v)
    b2=ap(2,1)*sin(u)+bp(2,1)*cos(u)+ap(2,2)*sin(v)+bp(2,2)*cos(v)
    f1=1+(a1-b1)**2+(a2-b2)**2
    f2=(u+3)**2+(v+1)**2
    y(i,:)=[-f1,-f2]
    end do
  end subroutine
  subroutine mop2(x,y)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::a
    integer::i,n
    n=size(x,2)
    a=1.0_dp/sqrt(real(n,dp))
    allocate(y(size(x,1),2))
    do i=1,size(x,1)
    y(i,1)=1-exp(-sum((4*x(i,:)-2-a)**2))
    y(i,2)=1-exp(-sum((4*x(i,:)-2+a)**2))
    end do
  end subroutine
  subroutine mop3(x,y)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::r2
    integer::i
    allocate(y(size(x,1),3))
    do i=1,size(x,1)
    r2=sum(x(i,1:2)**2)
    y(i,1)=0.5_dp*r2+sin(r2)
    y(i,2)=(3*x(i,1)-2*x(i,2)+4)**2/8+(x(i,1)-x(i,2)+1)**2/27+15
    y(i,3)=1/(r2+1)-1.1_dp*exp(-r2)
    end do
  end subroutine
  subroutine dtlz1(x,nobj,y)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::nobj
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::g,prodv
    integer::i,j,n,m
    n=size(x,2)
    m=nobj
    allocate(y(size(x,1),m))
    do i=1,size(x,1)
    g=100.0_dp*(real(n-m+1,dp)+sum((x(i,m:n)-0.5_dp)**2-cos(20*pi*(x(i,m:n)-0.5_dp))))
    do j=1,m
    prodv=1.0_dp
    if(m-j>=1)prodv=product(x(i,1:m-j))
    if(j>1)prodv=prodv*(1.0_dp-x(i,m-j+1))
    y(i,j)=0.5_dp*(1+g)*prodv
    end do
    end do
  end subroutine
  subroutine dtlz2(x,nobj,y)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::nobj
    real(dp),allocatable,intent(out)::y(:,:)
    call dtlz23(x,nobj,.false.,y)
  end subroutine
  subroutine dtlz3(x,nobj,y)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::nobj
    real(dp),allocatable,intent(out)::y(:,:)
    call dtlz23(x,nobj,.true.,y)
  end subroutine
  subroutine dtlz23(x,nobj,is3,y)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::nobj
    logical,intent(in)::is3
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::g,prodv
    integer::i,j,n,m
    n=size(x,2)
    m=nobj
    allocate(y(size(x,1),m))
    do i=1,size(x,1)
    if(is3)then
    g=100.0_dp*(real(n-m+1,dp)+sum((x(i,m:n)-0.5_dp)**2-cos(20*pi*(x(i,m:n)-0.5_dp))))
    else
    g=sum((x(i,m:n)-0.5_dp)**2)
    end if
    do j=1,m
    prodv=1.0_dp
    if(m-j>=1)prodv=product(cos(x(i,1:m-j)*pi/2))
    if(j>1)prodv=prodv*sin(x(i,m-j+1)*pi/2)
    y(i,j)=(1+g)*prodv
    end do
    end do
  end subroutine
  subroutine dtlz7(x,nobj,y)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::nobj
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::g,h
    integer::i,j,n,m
    n=size(x,2)
    m=nobj
    allocate(y(size(x,1),m))
    do i=1,size(x,1)
      g=1.0_dp
      do j=m,n
        g=g+9.0_dp*x(i,j)/real(j-m+1,dp)
      end do
      y(i,1:m-1)=x(i,1:m-1)
      h=real(m,dp)-sum(x(i,1:m-1)/(1+g)*(1+sin(3*pi*x(i,1:m-1))))
      y(i,m)=(1+g)*h
    end do
  end subroutine
  subroutine oka1(x,y)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::y(:,:)
    real(dp)::u,v,x1p,x2p,c,s
    integer::i
    c=cos(pi/12)
    s=sin(pi/12)
    allocate(y(size(x,1),2))
    do i=1,size(x,1)
    u=x(i,1)*2*pi*c+6*s
    v=x(i,2)*(6*c+2*pi*s)-2*pi*s
    x1p=c*u-s*v
    x2p=s*u+c*v
    y(i,1)=x1p
    y(i,2)=sqrt(2*pi)-sqrt(abs(x1p))+abs(x2p-3*cos(x1p)-3)**(1.0_dp/3.0_dp)
    end do
  end subroutine
end module gpareto_benchmarks
