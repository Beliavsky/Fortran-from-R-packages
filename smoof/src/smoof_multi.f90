! SPDX-License-Identifier: BSD-2-Clause
module smoof_multi
  use smoof_kinds, only : dp, pi
  implicit none
  private
  public :: dtlz1, dtlz2, dtlz3, dtlz4, dtlz5, dtlz6, dtlz7
  public :: zdt1, zdt2, zdt3, zdt4, zdt6
  public :: mop1, mop2, mop3, mop4, mop5, mop6, mop7
  public :: bk1, viennet, kursawe, dent, bi_sphere
  public :: wfg1, wfg2, wfg3, wfg4, wfg5, wfg6, wfg7, wfg8, wfg9

contains

  pure subroutine dtlz1(x, m, f)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: m
    real(dp), intent(out) :: f(m)
    integer :: i, n, k, j
    real(dp) :: g, prod_x
    n=size(x); k=n-m+1
    g=0.0_dp
    do j=m,n
      g=g+(x(j)-0.5_dp)**2-cos(20.0_dp*pi*(x(j)-0.5_dp))
    end do
    g=100.0_dp*(real(k,dp)+g)
    f=0.5_dp*(1.0_dp+g)
    prod_x=1.0_dp
    do i=m,2,-1
      f(i)=f(i)*prod_x*(1.0_dp-x(m-i+1))
      prod_x=prod_x*x(m-i+1)
    end do
    f(1)=f(1)*prod_x
  end subroutine dtlz1

  pure subroutine dtlz2(x, m, f)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: m
    real(dp), intent(out) :: f(m)
    integer :: i, j, n
    real(dp) :: g, prod_x
    n=size(x); g=0.0_dp
    do j=m,n; g=g+(x(j)-0.5_dp)**2; end do
    f=1.0_dp+g; prod_x=1.0_dp
    do i=m,2,-1
      f(i)=f(i)*prod_x*sin(x(m-i+1)*pi*0.5_dp)
      prod_x=prod_x*cos(x(m-i+1)*pi*0.5_dp)
    end do
    f(1)=f(1)*prod_x
  end subroutine dtlz2

  pure subroutine dtlz3(x, m, f)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: m
    real(dp), intent(out) :: f(m)
    integer :: i,j,n,k
    real(dp) :: g,prod_x
    n=size(x); k=n-m+1; g=0.0_dp
    do j=m,n
      g=g+(x(j)-0.5_dp)**2-cos(20.0_dp*pi*(x(j)-0.5_dp))
    end do
    g=100.0_dp*(real(k,dp)+g); f=1.0_dp+g; prod_x=1.0_dp
    do i=m,2,-1
      f(i)=f(i)*prod_x*sin(x(m-i+1)*pi*0.5_dp)
      prod_x=prod_x*cos(x(m-i+1)*pi*0.5_dp)
    end do
    f(1)=f(1)*prod_x
  end subroutine dtlz3

  pure subroutine dtlz4(x, m, f, alpha)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: m
    real(dp), intent(out) :: f(m)
    real(dp), intent(in), optional :: alpha
    integer :: i,j,n
    real(dp) :: g,prod_x,a
    real(dp) :: y(max(1,m-1))
    a=100.0_dp; if(present(alpha)) a=alpha
    n=size(x); g=0.0_dp
    do j=m,n; g=g+(x(j)-0.5_dp)**2; end do
    do j=1,m-1; y(j)=x(j)**a; end do
    f=1.0_dp+g; prod_x=1.0_dp
    do i=m,2,-1
      f(i)=f(i)*prod_x*sin(y(m-i+1)*pi*0.5_dp)
      prod_x=prod_x*cos(y(m-i+1)*pi*0.5_dp)
    end do
    f(1)=f(1)*prod_x
  end subroutine dtlz4

  pure subroutine dtlz5(x, m, f)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: m
    real(dp), intent(out) :: f(m)
    integer :: i,j,n
    real(dp) :: g,t,prod_x
    real(dp) :: theta(max(1,m-1))
    n=size(x); g=0.0_dp
    do j=m,n; g=g+(x(j)-0.5_dp)**2; end do
    theta(1)=x(1)*pi/2.0_dp
    t=pi/(4.0_dp*(1.0_dp+g))
    do j=2,m-1; theta(j)=t*(1.0_dp+2.0_dp*g*x(j)); end do
    f=1.0_dp+g; prod_x=1.0_dp
    do i=m,2,-1
      f(i)=f(i)*prod_x*sin(theta(m-i+1))
      prod_x=prod_x*cos(theta(m-i+1))
    end do
    f(1)=f(1)*prod_x
  end subroutine dtlz5

  pure subroutine dtlz6(x, m, f)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: m
    real(dp), intent(out) :: f(m)
    integer :: i,j,n
    real(dp) :: g,t,prod_x
    real(dp) :: theta(max(1,m-1))
    n=size(x); g=0.0_dp
    do j=m,n; g=g+x(j)**0.1_dp; end do
    theta(1)=x(1)*pi/2.0_dp
    t=pi/(4.0_dp*(1.0_dp+g))
    do j=2,m-1; theta(j)=t*(1.0_dp+2.0_dp*g*x(j)); end do
    f=1.0_dp+g; prod_x=1.0_dp
    do i=m,2,-1
      f(i)=f(i)*prod_x*sin(theta(m-i+1))
      prod_x=prod_x*cos(theta(m-i+1))
    end do
    f(1)=f(1)*prod_x
  end subroutine dtlz6

  pure subroutine dtlz7(x, m, f)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: m
    real(dp), intent(out) :: f(m)
    integer :: n,k,j
    real(dp) :: g,h
    n=size(x); k=n-m+1
    f(1:m-1)=x(1:m-1)
    g=1.0_dp+9.0_dp*sum(x(m:n))/real(k,dp)
    h=real(m,dp)
    do j=1,m-1
      h=h-f(j)*(1.0_dp+sin(3.0_dp*pi*f(j)))/(1.0_dp+g)
    end do
    f(m)=(1.0_dp+g)*h
  end subroutine dtlz7

  pure subroutine zdt1(x,f)
    real(dp), intent(in) :: x(:); real(dp), intent(out) :: f(2)
    real(dp) :: g
    g=1.0_dp+9.0_dp*sum(x(2:))/real(size(x)-1,dp)
    f(1)=x(1); f(2)=g*(1.0_dp-sqrt(f(1)/g))
  end subroutine zdt1

  pure subroutine zdt2(x,f)
    real(dp), intent(in) :: x(:); real(dp), intent(out) :: f(2)
    real(dp) :: g
    g=1.0_dp+9.0_dp*sum(x(2:))/real(size(x)-1,dp)
    f(1)=x(1); f(2)=g*(1.0_dp-(f(1)/g)**2)
  end subroutine zdt2

  pure subroutine zdt3(x,f)
    real(dp), intent(in) :: x(:); real(dp), intent(out) :: f(2)
    real(dp) :: g,r
    g=1.0_dp+9.0_dp*sum(x(2:))/real(size(x)-1,dp)
    f(1)=x(1); r=f(1)/g
    f(2)=g*(1.0_dp-sqrt(r)-r*sin(10.0_dp*pi*f(1)))
  end subroutine zdt3

  pure subroutine zdt4(x,f)
    real(dp), intent(in) :: x(:); real(dp), intent(out) :: f(2)
    real(dp) :: g
    g=1.0_dp+10.0_dp*real(size(x)-1,dp)+sum(x(2:)**2-10.0_dp*cos(4.0_dp*pi*x(2:)))
    f(1)=x(1); f(2)=g*(1.0_dp-sqrt(f(1)/g))
  end subroutine zdt4

  pure subroutine zdt6(x,f)
    real(dp), intent(in) :: x(:); real(dp), intent(out) :: f(2)
    real(dp) :: g
    f(1)=1.0_dp-exp(-4.0_dp*x(1))*sin(6.0_dp*pi*x(1))**6
    g=1.0_dp+9.0_dp*(sum(x(2:))/real(size(x)-1,dp))**0.25_dp
    f(2)=g*(1.0_dp-(f(1)/g)**2)
  end subroutine zdt6

  pure subroutine mop1(x,f)
    real(dp), intent(in)::x(:); real(dp),intent(out)::f(2)
    f=[x(1)**2,(x(1)-2.0_dp)**2]
  end subroutine mop1

  pure subroutine mop2(x,f)
    real(dp), intent(in)::x(:); real(dp),intent(out)::f(2)
    real(dp) :: r
    r=1.0_dp/sqrt(real(size(x),dp))
    f(1)=1.0_dp-exp(-sum((x-r)**2))
    f(2)=1.0_dp-exp(-sum((x+r)**2))
  end subroutine mop2

  pure subroutine mop3(x,f)
    real(dp), intent(in)::x(:); real(dp),intent(out)::f(2)
    real(dp)::a1,a2,b1,b2
    a1=0.5_dp*sin(1.0_dp)-2.0_dp*cos(1.0_dp)+sin(2.0_dp)-1.5_dp*cos(2.0_dp)
    a2=1.5_dp*sin(1.0_dp)-cos(1.0_dp)+2.0_dp*sin(2.0_dp)-0.5_dp*cos(2.0_dp)
    b1=0.5_dp*sin(x(1))-2.0_dp*cos(x(1))+sin(x(2))-1.5_dp*cos(x(2))
    b2=1.5_dp*sin(x(1))-cos(x(1))+2.0_dp*sin(x(2))-0.5_dp*cos(x(2))
    f(1)=-1.0_dp-(a1-b1)**2-(a2-b2)**2
    f(2)=-(x(1)+3.0_dp)**2-(x(2)+1.0_dp)**2
  end subroutine mop3

  pure subroutine mop4(x,f)
    real(dp), intent(in)::x(:); real(dp),intent(out)::f(2)
    integer::i
    f=0.0_dp
    do i=1,size(x)-1
      f(1)=f(1)-10.0_dp*exp(-0.2_dp*sqrt(x(i)**2+x(i+1)**2))
    end do
    f(2)=sum(abs(x)**0.8_dp+5.0_dp*sin(x**3))
  end subroutine mop4

  pure subroutine mop5(x,f)
    real(dp), intent(in)::x(:); real(dp),intent(out)::f(3)
    real(dp)::s
    s=x(1)**2+x(2)**2
    f(1)=0.5_dp*s+sin(s)
    f(2)=(3.0_dp*x(1)-2.0_dp*x(2)+4.0_dp)**2/8.0_dp &
      +(x(1)-x(2)+1.0_dp)**2/27.0_dp+15.0_dp
    f(3)=1.0_dp/(s+1.0_dp)-1.1_dp*exp(-s)
  end subroutine mop5

  pure subroutine mop6(x,f)
    real(dp), intent(in)::x(:); real(dp),intent(out)::f(2)
    real(dp)::a,b
    f(1)=x(1); a=1.0_dp+10.0_dp*x(2); b=x(1)/a
    f(2)=a*(1.0_dp-b*b-b*sin(8.0_dp*pi*x(1)))
  end subroutine mop6

  pure subroutine mop7(x,f)
    real(dp), intent(in)::x(:); real(dp),intent(out)::f(3)
    f(1)=(x(1)-2.0_dp)**2/2.0_dp+(x(2)+1.0_dp)**2/13.0_dp+3.0_dp
    f(2)=(x(1)+x(2)-3.0_dp)**2/36.0_dp+(-x(1)+x(2)+2.0_dp)**2/8.0_dp-17.0_dp
    f(3)=(x(1)+2.0_dp*x(2)-1.0_dp)**2/175.0_dp+(-x(1)+2.0_dp*x(2))**2/17.0_dp-13.0_dp
  end subroutine mop7

  pure subroutine bk1(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(2)
    f(1)=x(1)+x(2); f(2)=(x(1)-5.0_dp)**2+(x(2)-5.0_dp)**2
  end subroutine bk1

  pure subroutine viennet(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(3)
    real(dp)::s
    s=x(1)**2+x(2)**2
    f(1)=0.5_dp*s+sin(s)
    f(2)=(3.0_dp*x(1)-2.0_dp*x(2)+4.0_dp)**2/8.0_dp &
      +(x(1)-x(2)+1.0_dp)**2/27.0_dp+15.0_dp
    f(3)=1.0_dp/(s+1.0_dp)-1.1_dp*exp(-s)
  end subroutine viennet

  pure subroutine kursawe(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(2)
    integer::i
    f=0.0_dp
    do i=1,size(x)-1
      f(1)=f(1)-10.0_dp*exp(-0.2_dp*sqrt(x(i)**2+x(i+1)**2))
    end do
    f(2)=sum(abs(x)**0.8_dp+5.0_dp*sin(x)**3)
  end subroutine kursawe

  pure subroutine dent(x,f)
    real(dp),intent(in)::x(:);real(dp),intent(out)::f(2)
    real(dp)::d,a,b
    d=0.85_dp*exp(-(x(1)-x(2))**2)
    a=sqrt(1.0_dp+(x(1)+x(2))**2)
    b=sqrt(1.0_dp+(x(1)-x(2))**2)
    f(1)=0.5_dp*(a+b+x(1)-x(2))+d
    f(2)=0.5_dp*(a+b-x(1)+x(2))+d
  end subroutine dent

  pure subroutine bi_sphere(x,a,f)
    real(dp),intent(in)::x(:),a;real(dp),intent(out)::f(2)
    f(1)=sum(x*x); f(2)=sum((x-a)**2)
  end subroutine bi_sphere

  pure real(dp) function wfg_shift_linear(y,a) result(r)
    real(dp),intent(in)::y,a
    r=abs(y-a)/abs(real(floor(a-y),dp)+a)
  end function wfg_shift_linear

  pure real(dp) function wfg_bias_flat(y,a,b,c) result(r)
    real(dp),intent(in)::y,a,b,c
    r=-min(0.0_dp,real(floor(c-y),dp))*((1.0_dp-a)*(y-c)/(1.0_dp-c))
    r=r+a+min(0.0_dp,real(floor(y-b),dp))*(a*(b-y)/b)
  end function wfg_bias_flat

  pure real(dp) function wfg_bias_param(y,u,a,b,c) result(r)
    real(dp),intent(in)::y,u,a,b,c
    real(dp)::v
    v=a-(1.0_dp-2.0_dp*u)*abs(real(floor(0.5_dp-u),dp)+a)
    r=y**(b+(c-b)*v)
  end function wfg_bias_param

  pure real(dp) function wfg_shift_deceptive(y,a,b,c) result(r)
    real(dp),intent(in)::y,a,b,c
    r=real(floor(y-a+b),dp)*(1.0_dp-c+(a-b)/b)/(a-b)
    r=r+real(floor(a+b-y),dp)*(1.0_dp-c+(1.0_dp-a-b)/b)/(1.0_dp-a-b)
    r=1.0_dp+(abs(y-a)-b)*(r+1.0_dp/b)
  end function wfg_shift_deceptive

  pure real(dp) function wfg_shift_multimodal(y,a,b,c) result(r)
    real(dp),intent(in)::y,a,b,c
    real(dp)::q
    q=0.5_dp-abs(y-c)/(2.0_dp*(real(floor(c-y),dp)+c))
    r=(1.0_dp+cos((4.0_dp*a+2.0_dp)*pi*q) &
      +4.0_dp*b*(abs(y-c)/(2.0_dp*(real(floor(c-y),dp)+c)))**2)/(b+2.0_dp)
  end function wfg_shift_multimodal

  pure real(dp) function wfg_rsum(y,w) result(r)
    real(dp),intent(in)::y(:),w(:)
    r=sum(w*y)/sum(w)
  end function wfg_rsum

  pure real(dp) function wfg_rnonsep(y,a) result(r)
    real(dp),intent(in)::y(:);integer,intent(in)::a
    integer::j,k,n
    real(dp)::num,den
    n=size(y); num=0.0_dp
    do j=1,n
      num=num+y(j)
      do k=0,a-2
        num=num+abs(y(j)-y(mod(j+k,n)+1))
      end do
    end do
    den=(real(n,dp)/real(a,dp))*real(ceiling(real(a,dp)/2.0_dp),dp) &
      *(1.0_dp+2.0_dp*real(a,dp)-2.0_dp*real(ceiling(real(a,dp)/2.0_dp),dp))
    r=num/den
  end function wfg_rnonsep

  pure subroutine wfg_calc_x(t,a,x)
    real(dp),intent(in)::t(:),a(:);real(dp),intent(out)::x(size(t))
    integer::i,m
    m=size(t)
    do i=1,m-1
      x(i)=max(t(m),a(i))*(t(i)-0.5_dp)+0.5_dp
    end do
    x(m)=t(m)
  end subroutine wfg_calc_x

  pure real(dp) function wfg_shape_linear(x,m,obj) result(r)
    real(dp),intent(in)::x(:);integer,intent(in)::m,obj
    integer::i
    r=1.0_dp
    if(obj==1) then
      do i=1,m-1;r=r*x(i);end do
    else if(obj<m) then
      do i=1,m-obj;r=r*x(i);end do
      r=r*(1.0_dp-x(m-obj+1))
    else
      r=1.0_dp-x(1)
    end if
  end function wfg_shape_linear

  pure real(dp) function wfg_shape_convex(x,m,obj) result(r)
    real(dp),intent(in)::x(:);integer,intent(in)::m,obj
    integer::i
    r=1.0_dp
    if(obj==1) then
      do i=1,m-1;r=r*(1.0_dp-cos(x(i)*pi/2.0_dp));end do
    else if(obj<m) then
      do i=1,m-obj;r=r*(1.0_dp-cos(x(i)*pi/2.0_dp));end do
      r=r*(1.0_dp-sin(x(m-obj+1)*pi/2.0_dp))
    else
      r=1.0_dp-sin(x(1)*pi/2.0_dp)
    end if
  end function wfg_shape_convex

  pure real(dp) function wfg_shape_concave(x,m,obj) result(r)
    real(dp),intent(in)::x(:);integer,intent(in)::m,obj
    integer::i
    r=1.0_dp
    if(obj==1) then
      do i=1,m-1;r=r*sin(x(i)*pi/2.0_dp);end do
    else if(obj<m) then
      do i=1,m-obj;r=r*sin(x(i)*pi/2.0_dp);end do
      r=r*cos(x(m-obj+1)*pi/2.0_dp)
    else
      r=cos(x(1)*pi/2.0_dp)
    end if
  end function wfg_shape_concave

  pure real(dp) function wfg_shape_mixed(x,alpha,a) result(r)
    real(dp),intent(in)::x(:),alpha;integer,intent(in)::a
    r=(1.0_dp-x(1)-cos(2.0_dp*real(a,dp)*pi*x(1)+pi/2.0_dp) &
      /(2.0_dp*real(a,dp)*pi))**alpha
  end function wfg_shape_mixed

  pure real(dp) function wfg_shape_disc(x,alpha,beta,a) result(r)
    real(dp),intent(in)::x(:),alpha,beta;integer,intent(in)::a
    r=1.0_dp-x(1)**alpha*cos(real(a,dp)*x(1)**beta*pi)**2
  end function wfg_shape_disc

  pure subroutine wfg_finish(t,a,m,kind,f)
    real(dp),intent(in)::t(:),a(:);integer,intent(in)::m,kind
    real(dp),intent(out)::f(m)
    real(dp)::x(m),h
    integer::j
    call wfg_calc_x(t,a,x)
    do j=1,m
      select case(kind)
      case(1); h=wfg_shape_linear(x,m,j)
      case(2); h=wfg_shape_convex(x,m,j)
      case default; h=wfg_shape_concave(x,m,j)
      end select
      f(j)=x(m)+2.0_dp*real(j,dp)*h
    end do
  end subroutine wfg_finish

  pure subroutine wfg1(z,m,k,f)
    real(dp),intent(in)::z(:);integer,intent(in)::m,k;real(dp),intent(out)::f(m)
    real(dp)::y(size(z)),w(size(z)),t(m),a(max(1,m-1)),x(m),h
    integer::i,j,tmp,n
    n=size(z); a=1.0_dp
    do i=1,n;y(i)=z(i)/(2.0_dp*real(i,dp));w(i)=2.0_dp*real(i,dp);end do
    do i=k+1,n
      y(i)=wfg_shift_linear(y(i),0.35_dp)
      y(i)=wfg_bias_flat(y(i),0.8_dp,0.75_dp,0.85_dp)
    end do
    y=y**0.02_dp; tmp=k/(m-1)
    do i=1,m-1
      t(i)=wfg_rsum(y((i-1)*tmp+1:i*tmp),w((i-1)*tmp+1:i*tmp))
    end do
    t(m)=wfg_rsum(y(k+1:n),w(k+1:n))
    call wfg_calc_x(t,a,x)
    do j=1,m-1
      h=wfg_shape_convex(x,m,j);f(j)=x(m)+2.0_dp*real(j,dp)*h
    end do
    f(m)=x(m)+2.0_dp*real(m,dp)*wfg_shape_mixed(x,1.0_dp,5)
  end subroutine wfg1

  pure subroutine wfg2(z,m,k,f)
    real(dp),intent(in)::z(:);integer,intent(in)::m,k;real(dp),intent(out)::f(m)
    real(dp)::y(size(z)),t2(k+(size(z)-k)/2),t3(m),a(max(1,m-1)),x(m),h
    real(dp),allocatable::ones(:)
    integer::i,j,n,l2,tmp
    n=size(z);l2=(n-k)/2;a=1.0_dp
    do i=1,n;y(i)=z(i)/(2.0_dp*real(i,dp));end do
    do i=k+1,n;y(i)=wfg_shift_linear(y(i),0.35_dp);end do
    t2(1:k)=y(1:k)
    do i=1,l2;t2(k+i)=wfg_rnonsep(y(k+2*i-1:k+2*i),2);end do
    tmp=k/(m-1)
    allocate(ones(max(tmp,l2)));ones=1.0_dp
    do i=1,m-1;t3(i)=wfg_rsum(t2((i-1)*tmp+1:i*tmp),ones(1:tmp));end do
    t3(m)=wfg_rsum(t2(k+1:k+l2),ones(1:l2));deallocate(ones)
    call wfg_calc_x(t3,a,x)
    do j=1,m-1
      h=wfg_shape_convex(x,m,j);f(j)=x(m)+2.0_dp*real(j,dp)*h
    end do
    f(m)=x(m)+2.0_dp*real(m,dp)*wfg_shape_disc(x,1.0_dp,1.0_dp,5)
  end subroutine wfg2

  pure subroutine wfg3(z,m,k,f)
    real(dp),intent(in)::z(:);integer,intent(in)::m,k;real(dp),intent(out)::f(m)
    real(dp)::y(size(z)),t2(k+(size(z)-k)/2),t3(m),a(max(1,m-1))
    real(dp),allocatable::ones(:)
    integer::i,n,l2,tmp
    n=size(z);l2=(n-k)/2;a=0.0_dp;a(1)=1.0_dp
    do i=1,n;y(i)=z(i)/(2.0_dp*real(i,dp));end do
    do i=k+1,n;y(i)=wfg_shift_linear(y(i),0.35_dp);end do
    t2(1:k)=y(1:k)
    do i=1,l2;t2(k+i)=wfg_rnonsep(y(k+2*i-1:k+2*i),2);end do
    tmp=k/(m-1);allocate(ones(max(tmp,l2)));ones=1.0_dp
    do i=1,m-1;t3(i)=wfg_rsum(t2((i-1)*tmp+1:i*tmp),ones(1:tmp));end do
    t3(m)=wfg_rsum(t2(k+1:k+l2),ones(1:l2));deallocate(ones)
    call wfg_finish(t3,a,m,1,f)
  end subroutine wfg3

  pure subroutine wfg4(z,m,k,f)
    real(dp),intent(in)::z(:);integer,intent(in)::m,k;real(dp),intent(out)::f(m)
    real(dp)::y(size(z)),t(m),a(max(1,m-1))
    real(dp),allocatable::ones(:)
    integer::i,n,tmp
    n=size(z);a=1.0_dp
    do i=1,n
      y(i)=z(i)/(2.0_dp*real(i,dp));y(i)=wfg_shift_multimodal(y(i),30.0_dp,10.0_dp,0.35_dp)
    end do
    tmp=k/(m-1);allocate(ones(max(tmp,n-k)));ones=1.0_dp
    do i=1,m-1;t(i)=wfg_rsum(y((i-1)*tmp+1:i*tmp),ones(1:tmp));end do
    t(m)=wfg_rsum(y(k+1:n),ones(1:n-k));deallocate(ones)
    call wfg_finish(t,a,m,3,f)
  end subroutine wfg4

  pure subroutine wfg5(z,m,k,f)
    real(dp),intent(in)::z(:);integer,intent(in)::m,k;real(dp),intent(out)::f(m)
    real(dp)::y(size(z)),t(m),a(max(1,m-1))
    real(dp),allocatable::ones(:)
    integer::i,n,tmp
    n=size(z);a=1.0_dp
    do i=1,n
      y(i)=z(i)/(2.0_dp*real(i,dp));y(i)=wfg_shift_deceptive(y(i),0.35_dp,0.001_dp,0.05_dp)
    end do
    tmp=k/(m-1);allocate(ones(max(tmp,n-k)));ones=1.0_dp
    do i=1,m-1;t(i)=wfg_rsum(y((i-1)*tmp+1:i*tmp),ones(1:tmp));end do
    t(m)=wfg_rsum(y(k+1:n),ones(1:n-k));deallocate(ones)
    call wfg_finish(t,a,m,3,f)
  end subroutine wfg5

  pure subroutine wfg6(z,m,k,f)
    real(dp),intent(in)::z(:);integer,intent(in)::m,k;real(dp),intent(out)::f(m)
    real(dp)::y(size(z)),t(m),a(max(1,m-1))
    integer::i,n,l,tmp
    n=size(z);l=n-k;a=1.0_dp
    do i=1,n;y(i)=z(i)/(2.0_dp*real(i,dp));end do
    do i=k+1,n;y(i)=wfg_shift_linear(y(i),0.35_dp);end do
    tmp=k/(m-1)
    do i=1,m-1;t(i)=wfg_rnonsep(y((i-1)*tmp+1:i*tmp),tmp);end do
    t(m)=wfg_rnonsep(y(k+1:n),l)
    call wfg_finish(t,a,m,3,f)
  end subroutine wfg6

  pure subroutine wfg7(z,m,k,f)
    real(dp),intent(in)::z(:);integer,intent(in)::m,k;real(dp),intent(out)::f(m)
    real(dp)::y(size(z)),t(m),a(max(1,m-1)),u
    real(dp),allocatable::ones(:)
    integer::i,n,l,tmp
    n=size(z);l=n-k;a=1.0_dp
    do i=1,n;y(i)=z(i)/(2.0_dp*real(i,dp));end do
    allocate(ones(n));ones=1.0_dp
    do i=1,k
      u=wfg_rsum(y(i+1:n),ones(1:n-i))
      y(i)=wfg_bias_param(y(i),u,0.98_dp/49.98_dp,0.02_dp,50.0_dp)
    end do
    do i=k+1,n;y(i)=wfg_shift_linear(y(i),0.35_dp);end do
    tmp=k/(m-1)
    do i=1,m-1;t(i)=wfg_rsum(y((i-1)*tmp+1:i*tmp),ones(1:tmp));end do
    t(m)=wfg_rsum(y(k+1:n),ones(1:l));deallocate(ones)
    call wfg_finish(t,a,m,3,f)
  end subroutine wfg7

  pure subroutine wfg8(z,m,k,f)
    real(dp),intent(in)::z(:);integer,intent(in)::m,k;real(dp),intent(out)::f(m)
    real(dp)::y(size(z)),t(m),a(max(1,m-1)),u
    real(dp),allocatable::ones(:)
    integer::i,n,l,tmp
    n=size(z);l=n-k;a=1.0_dp
    do i=1,n;y(i)=z(i)/(2.0_dp*real(i,dp));end do
    allocate(ones(n));ones=1.0_dp
    do i=n,k+1,-1
      u=wfg_rsum(y(1:i-1),ones(1:i-1))
      y(i)=wfg_bias_param(y(i),u,0.98_dp/49.98_dp,0.02_dp,50.0_dp)
    end do
    do i=k+1,n;y(i)=wfg_shift_linear(y(i),0.35_dp);end do
    tmp=k/(m-1)
    do i=1,m-1;t(i)=wfg_rsum(y((i-1)*tmp+1:i*tmp),ones(1:tmp));end do
    t(m)=wfg_rsum(y(k+1:n),ones(1:l));deallocate(ones)
    call wfg_finish(t,a,m,3,f)
  end subroutine wfg8

  pure subroutine wfg9(z,m,k,f)
    real(dp),intent(in)::z(:);integer,intent(in)::m,k;real(dp),intent(out)::f(m)
    real(dp)::y(size(z)),t(m),a(max(1,m-1)),u
    real(dp),allocatable::ones(:)
    integer::i,n,l,tmp
    n=size(z);l=n-k;a=1.0_dp
    do i=1,n;y(i)=z(i)/(2.0_dp*real(i,dp));end do
    allocate(ones(n));ones=1.0_dp
    do i=1,n-1
      u=wfg_rsum(y(i+1:n),ones(1:n-i))
      y(i)=wfg_bias_param(y(i),u,0.98_dp/49.98_dp,0.02_dp,50.0_dp)
    end do
    do i=1,k;y(i)=wfg_shift_deceptive(y(i),0.35_dp,0.001_dp,0.05_dp);end do
    do i=k+1,n;y(i)=wfg_shift_multimodal(y(i),30.0_dp,95.0_dp,0.35_dp);end do
    tmp=k/(m-1)
    do i=1,m-1;t(i)=wfg_rnonsep(y((i-1)*tmp+1:i*tmp),tmp);end do
    t(m)=wfg_rnonsep(y(k+1:n),l);deallocate(ones)
    call wfg_finish(t,a,m,3,f)
  end subroutine wfg9

end module smoof_multi
