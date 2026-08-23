module kernsmooth_mod
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter :: pi = acos(-1.0_dp)
  public :: linbin, linbin2d, rlbin, bkde, bkde2d, bkfe, locpoly
  public :: dpih, dpik, dpill, sdiag, sstdiag, mean_dp, sd_dp, quantile_dp
contains

  pure real(dp) function mean_dp(x) result(v)
    real(dp), intent(in) :: x(:)
    v = sum(x)/real(size(x),dp)
  end function

  pure real(dp) function sd_dp(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    m = mean_dp(x)
    if (size(x) > 1) then
      v = sqrt(sum((x-m)**2)/real(size(x)-1,dp))
    else
      v = 0.0_dp
    end if
  end function

  function quantile_dp(x, p) result(q)
    real(dp), intent(in) :: x(:), p
    real(dp) :: q, h, frac
    real(dp), allocatable :: z(:)
    integer :: i, j, n, k
    n = size(x); allocate(z(n)); z=x
    do i=2,n
      h=z(i); j=i-1
      do while (j>=1 .and. z(j)>h)
        z(j+1)=z(j); j=j-1
      end do
      z(j+1)=h
    end do
    if (p<=0.0_dp) then; q=z(1); return; end if
    if (p>=1.0_dp) then; q=z(n); return; end if
    h=1.0_dp+(real(n-1,dp))*p; k=int(floor(h)); frac=h-real(k,dp)
    q=(1.0_dp-frac)*z(k)+frac*z(min(k+1,n))
  end function

  subroutine linbin(x, a, b, counts, truncate)
    real(dp), intent(in) :: x(:), a, b
    real(dp), intent(out) :: counts(:)
    logical, intent(in), optional :: truncate
    logical :: tr
    real(dp) :: u, w
    integer :: i, j, m
    m=size(counts); counts=0.0_dp; tr=.true.; if(present(truncate)) tr=truncate
    if (m<2 .or. b<=a) return
    do i=1,size(x)
      u=(x(i)-a)*real(m-1,dp)/(b-a)+1.0_dp
      if (u<1.0_dp) then
        if (.not.tr) counts(1)=counts(1)+1.0_dp
      else if (u>real(m,dp)) then
        if (.not.tr) counts(m)=counts(m)+1.0_dp
      else
        j=min(int(floor(u)),m-1); w=u-real(j,dp)
        if (j<1) then
          counts(1)=counts(1)+1.0_dp
        else
          counts(j)=counts(j)+(1.0_dp-w)
          counts(j+1)=counts(j+1)+w
        end if
      end if
    end do
  end subroutine

  subroutine rlbin(x, y, a, b, xcounts, ycounts, truncate)
    real(dp), intent(in) :: x(:), y(:), a, b
    real(dp), intent(out) :: xcounts(:), ycounts(:)
    logical, intent(in), optional :: truncate
    logical :: tr
    real(dp) :: u,w
    integer :: i,j,m
    m=size(xcounts); xcounts=0.0_dp; ycounts=0.0_dp; tr=.true.; if(present(truncate)) tr=truncate
    do i=1,min(size(x),size(y))
      u=(x(i)-a)*real(m-1,dp)/(b-a)+1.0_dp
      if (u<1.0_dp) then
        if(.not.tr) then; xcounts(1)=xcounts(1)+1; ycounts(1)=ycounts(1)+y(i); end if
      else if (u>real(m,dp)) then
        if(.not.tr) then; xcounts(m)=xcounts(m)+1; ycounts(m)=ycounts(m)+y(i); end if
      else
        j=min(max(int(floor(u)),1),m-1); w=u-real(j,dp)
        xcounts(j)=xcounts(j)+1-w; xcounts(j+1)=xcounts(j+1)+w
        ycounts(j)=ycounts(j)+(1-w)*y(i); ycounts(j+1)=ycounts(j+1)+w*y(i)
      end if
    end do
  end subroutine

  subroutine linbin2d(x, a, b, counts)
    real(dp), intent(in) :: x(:,:), a(2), b(2)
    real(dp), intent(out) :: counts(:,:)
    real(dp) :: u1,u2,w1,w2
    integer :: i,j,k,m1,m2
    m1=size(counts,1);m2=size(counts,2);counts=0
    do i=1,size(x,1)
      u1=(x(i,1)-a(1))*real(m1-1,dp)/(b(1)-a(1))+1
      u2=(x(i,2)-a(2))*real(m2-1,dp)/(b(2)-a(2))+1
      if(u1<1.or.u1>m1.or.u2<1.or.u2>m2) cycle
      j=min(max(int(floor(u1)),1),m1-1); k=min(max(int(floor(u2)),1),m2-1)
      w1=u1-j;w2=u2-k
      counts(j,k)=counts(j,k)+(1-w1)*(1-w2)
      counts(j+1,k)=counts(j+1,k)+w1*(1-w2)
      counts(j,k+1)=counts(j,k+1)+(1-w1)*w2
      counts(j+1,k+1)=counts(j+1,k+1)+w1*w2
    end do
  end subroutine

  pure real(dp) function kernel_value(u, kernel) result(k)
    real(dp), intent(in) :: u
    character(len=*), intent(in) :: kernel
    real(dp) :: z
    z=abs(u)
    select case(trim(kernel))
    case('normal'); k=exp(-0.5_dp*u*u)/sqrt(2*pi)
    case('box'); if(z<=1) then;k=0.5_dp;else;k=0;end if
    case('epanech'); if(z<=1) then;k=0.75_dp*(1-u*u);else;k=0;end if
    case('biweight'); if(z<=1) then;k=15.0_dp/16*(1-u*u)**2;else;k=0;end if
    case('triweight'); if(z<=1) then;k=35.0_dp/32*(1-u*u)**3;else;k=0;end if
    case default; k=exp(-0.5_dp*u*u)/sqrt(2*pi)
    end select
  end function

  subroutine bkde(x, bandwidth, grid, dens, kernel, range_x, canonical, truncate)
    real(dp), intent(in) :: x(:), bandwidth
    real(dp), intent(out) :: grid(:), dens(:)
    character(len=*), intent(in), optional :: kernel
    real(dp), intent(in), optional :: range_x(2)
    logical, intent(in), optional :: canonical, truncate
    character(len=16) :: ker
    real(dp) :: a,b,h,del0,dx,tot
    real(dp), allocatable :: counts(:)
    integer :: i,j,m,n
    ker='normal';if(present(kernel))ker=kernel;m=size(grid);n=size(x)
    del0=1
    select case(trim(ker))
    case('normal');del0=(1/(4*pi))**0.1_dp
    case('box');del0=(9.0_dp/2)**0.2_dp
    case('epanech');del0=15.0_dp**0.2_dp
    case('biweight');del0=35.0_dp**0.2_dp
    case('triweight');del0=(9450.0_dp/143)**0.2_dp
    end select
    h=bandwidth;if(present(canonical)) then;if(canonical)h=del0*bandwidth;end if
    if(present(range_x))then;a=range_x(1);b=range_x(2);else;a=minval(x)-4*h;b=maxval(x)+4*h;end if
    dx=(b-a)/real(m-1,dp); do i=1,m;grid(i)=a+(i-1)*dx;end do
    allocate(counts(m)); call linbin(x,a,b,counts,truncate)
    dens=0
    do i=1,m
      do j=1,m
        dens(i)=dens(i)+counts(j)*kernel_value((grid(i)-grid(j))/h,ker)/(real(n,dp)*h)
      end do
    end do
    tot=sum(dens)*dx
    if(tot>0) dens=dens/tot
  end subroutine

  subroutine bkde2d(x, bandwidth, grid1, grid2, fhat, range_x)
    real(dp), intent(in) :: x(:,:), bandwidth(2)
    real(dp), intent(out) :: grid1(:),grid2(:),fhat(:,:)
    real(dp), intent(in), optional :: range_x(2,2)
    real(dp) :: a(2),b(2),d1,d2,h1,h2,tot
    real(dp),allocatable :: c(:,:)
    integer :: i,j,k,l,m1,m2,n
    m1=size(grid1);m2=size(grid2);n=size(x,1);h1=bandwidth(1);h2=bandwidth(2)
    if(present(range_x))then;a=range_x(:,1);b=range_x(:,2)
    else;a=[minval(x(:,1))-1.5_dp*h1,minval(x(:,2))-1.5_dp*h2];b=[maxval(x(:,1))+1.5_dp*h1,maxval(x(:,2))+1.5_dp*h2];end if
    d1=(b(1)-a(1))/real(m1-1,dp);d2=(b(2)-a(2))/real(m2-1,dp)
    do i=1,m1;grid1(i)=a(1)+(i-1)*d1;end do;do j=1,m2;grid2(j)=a(2)+(j-1)*d2;end do
    allocate(c(m1,m2));call linbin2d(x,a,b,c);fhat=0
    do i=1,m1;do j=1,m2;do k=1,m1;do l=1,m2
      fhat(i,j)=fhat(i,j)+c(k,l)*kernel_value((grid1(i)-grid1(k))/h1,'normal')* &
        kernel_value((grid2(j)-grid2(l))/h2,'normal')/(real(n,dp)*h1*h2)
    end do;end do;end do;end do
    where(fhat<0)fhat=0;tot=sum(fhat)*d1*d2;if(tot>0)fhat=fhat/tot
  end subroutine

  pure real(dp) function hermite_prob(n,x) result(h)
    integer,intent(in)::n;real(dp),intent(in)::x
    real(dp)::h0,h1,hn;integer::i
    if(n==0)then;h=1;return;else if(n==1)then;h=x;return;end if
    h0=1;h1=x
    do i=2,n;hn=x*h1-real(i-1,dp)*h0;h0=h1;h1=hn;end do
    h=h1
  end function

  function bkfe(x, drv, bandwidth, range_x, binned) result(v)
    real(dp),intent(in)::x(:),bandwidth,range_x(2);integer,intent(in)::drv
    logical,intent(in),optional::binned
    real(dp)::v,a,b,delta,arg,k,n
    real(dp),allocatable::c(:),g(:);integer::i,j,m
    logical::bin
    m=size(x);a=range_x(1);b=range_x(2);bin=.false.;if(present(binned))bin=binned
    if(bin)then;allocate(c(m));c=x;else;allocate(c(401));m=401;call linbin(x,a,b,c);end if
    allocate(g(m));delta=(b-a)/real(m-1,dp);do i=1,m;g(i)=a+(i-1)*delta;end do
    n=sum(c);v=0
    do i=1,m;do j=1,m
      arg=(g(i)-g(j))/bandwidth
      k=hermite_prob(drv,arg)*exp(-0.5_dp*arg*arg)/sqrt(2*pi)/bandwidth**(drv+1)
      v=v+c(i)*c(j)*k
    end do;end do
    v=v/(n*n)
  end function

  subroutine solve_linear(a,b,x,ok)
    real(dp),intent(in)::a(:,:),b(:);real(dp),intent(out)::x(:);logical,intent(out)::ok
    real(dp),allocatable::aa(:,:),bb(:);real(dp)::fac,tmp;integer::n,i,j,k,p
    n=size(b);allocate(aa(n,n),bb(n));aa=a;bb=b;ok=.true.
    do k=1,n-1;p=k;do i=k+1,n;if(abs(aa(i,k))>abs(aa(p,k)))p=i;end do
      if(abs(aa(p,k))<1e-14_dp)then;ok=.false.;x=0;return;end if
      if(p/=k)then;do j=k,n;tmp=aa(k,j);aa(k,j)=aa(p,j);aa(p,j)=tmp;end do;tmp=bb(k);bb(k)=bb(p);bb(p)=tmp;end if
      do i=k+1,n;fac=aa(i,k)/aa(k,k);aa(i,k:n)=aa(i,k:n)-fac*aa(k,k:n);bb(i)=bb(i)-fac*bb(k);end do
    end do
    if(abs(aa(n,n))<1e-14_dp)then;ok=.false.;x=0;return;end if
    x(n)=bb(n)/aa(n,n);do i=n-1,1,-1;x(i)=(bb(i)-dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i);end do
  end subroutine

  subroutine locpoly(x,y,bandwidth,grid,fit,drv,degree,range_x)
    real(dp),intent(in)::x(:),y(:),bandwidth;real(dp),intent(out)::grid(:),fit(:)
    integer,intent(in),optional::drv,degree;real(dp),intent(in),optional::range_x(2)
    integer::d,p,m,i,j,r,s,n;real(dp)::a,b,dx,u,w,fac
    real(dp),allocatable::sm(:,:),tv(:),coef(:);logical::ok
    d=0;if(present(drv))d=drv;p=d+1;if(present(degree))p=degree+1;m=size(grid);n=size(x)
    if(present(range_x))then;a=range_x(1);b=range_x(2);else;a=minval(x);b=maxval(x);end if
    dx=(b-a)/real(m-1,dp);do i=1,m;grid(i)=a+(i-1)*dx;end do
    allocate(sm(p,p),tv(p),coef(p));fit=0
    do i=1,m;sm=0;tv=0
      do j=1,n;u=(x(j)-grid(i))/bandwidth;w=kernel_value(u,'normal')
        do r=0,p-1;tv(r+1)=tv(r+1)+w*(x(j)-grid(i))**r*y(j);do s=0,p-1;sm(r+1,s+1)=sm(r+1,s+1)+w*(x(j)-grid(i))**(r+s);end do;end do
      end do
      call solve_linear(sm,tv,coef,ok);if(ok)then;fac=1;do j=2,d;fac=fac*j;end do;if(d==0)fac=1;fit(i)=fac*coef(d+1);end if
    end do
  end subroutine

  function scale_est(x,mode) result(s)
    real(dp),intent(in)::x(:);character(len=*),intent(in)::mode;real(dp)::s,iq
    iq=(quantile_dp(x,.75_dp)-quantile_dp(x,.25_dp))/1.349_dp
    select case(trim(mode));case('stdev');s=sd_dp(x);case('iqr');s=iq;case default;s=min(iq,sd_dp(x));end select
  end function

  function dpih(x, level, scalest) result(h)
    real(dp),intent(in)::x(:);integer,intent(in),optional::level;character(len=*),intent(in),optional::scalest
    real(dp)::h,s,alpha,psi2,psi4,a,b;real(dp),allocatable::sx(:),c(:);integer::lev,n,m
    character(len=16)::sc
    lev=2;if(present(level))lev=level;sc='minim';if(present(scalest))sc=scalest;n=size(x);s=scale_est(x,sc)
    allocate(sx(n));sx=(x-mean_dp(x))/s;a=minval(sx);b=maxval(sx);m=401;allocate(c(m));call linbin(sx,a,b,c)
    if(lev==0)then;h=(24*sqrt(pi)/n)**(1.0_dp/3);else
      alpha=(2.0_dp/(5*n))**(1.0_dp/7)*sqrt(2.0_dp);psi4=bkfe(c,4,alpha,[a,b],.true.)
      alpha=(sqrt(2.0_dp/pi)/(psi4*n))**(1.0_dp/5);psi2=bkfe(c,2,alpha,[a,b],.true.)
      h=(6.0_dp/(-psi2*n))**(1.0_dp/3)
    end if
    h=s*h
  end function

  function dpik(x, level, kernel, canonical, scalest) result(h)
    real(dp),intent(in)::x(:);integer,intent(in),optional::level;character(len=*),intent(in),optional::kernel,scalest
    logical,intent(in),optional::canonical
    real(dp)::h,s,alpha,psi4,a,b,del0;real(dp),allocatable::sx(:),c(:);integer::lev,n,m
    character(len=16)::ker,sc
    lev=2;if(present(level))lev=level;ker='normal';if(present(kernel))ker=kernel;sc='minim';if(present(scalest))sc=scalest
    n = size(x)
    s = scale_est(x,sc)
    allocate(sx(n))
    sx = (x-mean_dp(x))/s
    a = minval(sx); b = maxval(sx); m = 401
    allocate(c(m))
    call linbin(sx,a,b,c)
    if (lev == 0) then
      psi4 = 3.0_dp/(8*sqrt(pi))
    else
      alpha = (2.0_dp*sqrt(2.0_dp)**9/(7*n))**(1.0_dp/9)
      psi4 = bkfe(c,6,alpha,[a,b],.true.)
      alpha = (-3*sqrt(2.0_dp/pi)/(psi4*n))**(1.0_dp/7)
      psi4 = bkfe(c,4,alpha,[a,b],.true.)
    end if
    select case(trim(ker))
    case('normal'); del0 = 1/(4*pi)**0.1_dp
    case('box'); del0 = (9.0_dp/2)**0.2_dp
    case('epanech'); del0 = 15.0_dp**0.2_dp
    case('biweight'); del0 = 35.0_dp**0.2_dp
    case('triweight'); del0 = (9450.0_dp/143)**0.2_dp
    case default; del0 = 1/(4*pi)**0.1_dp
    end select
    if(present(canonical))then;if(canonical)del0=1;end if
    h=s*del0*(1.0_dp/(psi4*n))**.2_dp
  end function

  subroutine local_linear_hat(x, bandwidth, grid, diagv, sqdiag)
    real(dp), intent(in) :: x(:), bandwidth
    real(dp), intent(out) :: grid(:), diagv(:), sqdiag(:)
    integer :: i, j, m, n
    real(dp) :: a, b, dx, u, w, s0, s1, s2, den, hij
    n = size(x); m = size(grid)
    a = minval(x); b = maxval(x)
    dx = (b-a)/real(m-1,dp)
    do i=1,m
      grid(i)=a+real(i-1,dp)*dx
      s0=0; s1=0; s2=0
      do j=1,n
        u=x(j)-grid(i); w=kernel_value(u/bandwidth,'normal')
        s0=s0+w; s1=s1+w*u; s2=s2+w*u*u
      end do
      den=s0*s2-s1*s1
      diagv(i)=0; sqdiag(i)=0
      if (abs(den)>1e-14_dp) then
        do j=1,n
          u=x(j)-grid(i); w=kernel_value(u/bandwidth,'normal')
          hij=w*(s2-s1*u)/den
          sqdiag(i)=sqdiag(i)+hij*hij
          if(abs(x(j)-grid(i))<=0.5_dp*dx) diagv(i)=diagv(i)+hij
        end do
      end if
    end do
  end subroutine

  subroutine sdiag(x, bandwidth, grid, values)
    real(dp), intent(in) :: x(:), bandwidth
    real(dp), intent(out) :: grid(:), values(:)
    real(dp), allocatable :: tmp(:)
    allocate(tmp(size(values)))
    call local_linear_hat(x,bandwidth,grid,values,tmp)
  end subroutine

  subroutine sstdiag(x, bandwidth, grid, values)
    real(dp), intent(in) :: x(:), bandwidth
    real(dp), intent(out) :: grid(:), values(:)
    real(dp), allocatable :: tmp(:)
    allocate(tmp(size(values)))
    call local_linear_hat(x,bandwidth,grid,tmp,values)
  end subroutine

  function dpill(x, y, trim) result(h)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: trim
    real(dp) :: h, tr, a, b, pilot, sig2, th22, dx
    real(dp), allocatable :: grid(:), fit(:), d2(:), resid(:)
    integer :: n, m, i, lo, hi
    n=min(size(x),size(y)); m=401
    tr=0.01_dp; if(present(trim)) tr=trim
    a=minval(x(1:n)); b=maxval(x(1:n))
    pilot=1.06_dp*sd_dp(x(1:n))*real(n,dp)**(-0.2_dp)
    if(pilot<=0) pilot=max((b-a)/20.0_dp,epsilon(1.0_dp))
    allocate(grid(m),fit(m),d2(m),resid(n))
    call locpoly(x(1:n),y(1:n),pilot,grid,fit,drv=0,degree=1,range_x=[a,b])
    do i=1,n
      resid(i)=y(i)-interp_linear(grid,fit,x(i))
    end do
    sig2=sum(resid*resid)/real(max(1,n-2),dp)
    call locpoly(x(1:n),y(1:n),pilot,grid,d2,drv=2,degree=3,range_x=[a,b])
    lo=max(1,1+int(tr*m)); hi=min(m,m-int(tr*m)); dx=(b-a)/real(m-1,dp)
    th22=sum(d2(lo:hi)**2)*dx/max(b-a,epsilon(1.0_dp))
    if(th22<=epsilon(1.0_dp)) then
      h=pilot
    else
      h=(sig2*(b-a)/(2*sqrt(pi)*th22*real(n,dp)))**0.2_dp
    end if
  end function

  pure real(dp) function interp_linear(xg, yg, x) result(v)
    real(dp), intent(in) :: xg(:), yg(:), x
    integer :: j, m
    real(dp) :: t
    m=size(xg)
    if(x<=xg(1)) then; v=yg(1); return; end if
    if(x>=xg(m)) then; v=yg(m); return; end if
    j=1
    do while(j<m-1 .and. x>xg(j+1)); j=j+1; end do
    t=(x-xg(j))/(xg(j+1)-xg(j))
    v=(1-t)*yg(j)+t*yg(j+1)
  end function

end module
