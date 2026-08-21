! Derived from locfit src/band.c, GPL-2-or-later.
module locfit_bandwidth
  use locfit_kinds, only : dp
  use locfit_constants
  use locfit_kernels, only : kernel_weight, kernel_derivative_over_u, &
    kernel_integral_moment, kernel_convolution, kernel_convolution1, &
    kernel_convolution4, kernel_convolution5, kernel_convolution6, kernel_wikk
  implicit none
  private
  public :: compsda, width_sj, kde_criterion, kde_select
  integer, parameter, public :: kde_aic=1, kde_ocv=2, kde_lscv=3
  integer, parameter, public :: kde_bcv=4, kde_sjpi=5, kde_gkk=6

contains

  pure subroutine sort_values(a)
    real(dp),intent(inout)::a(:)
    integer::i,j
    real(dp)::v
    do i=2,size(a)
      v=a(i);j=i-1
      do while(j>=1)
        if(a(j)<=v)exit
        a(j+1)=a(j);j=j-1
      end do
      a(j+1)=v
    end do
  end subroutine sort_values

  pure real(dp) function compsda(x,h) result(sd)
    real(dp),intent(in)::x(:),h
    real(dp)::ik,z
    integer::i,j,n
    n=size(x);sd=0.0_dp
    if(n<2 .or. h<=0.0_dp)return
    ik=kernel_integral_moment(1,ker=wgaus)
    do i=1,n
      do j=i,n
        z=(x(i)-x(j))/h
        sd=sd+merge(1.0_dp,2.0_dp,i==j)*kernel_convolution4(z,wgaus)/(ik*ik)
      end do
    end do
    sd=sd/(real(n*(n-1),dp)*h**5)
  end function compsda

  pure real(dp) function width_sj(x,lambda) result(c)
    real(dp),intent(in)::x(:),lambda
    real(dp)::ik,a,b,td,sa,z,c1,c2,c3
    integer::i,j,n
    integer::pow2(1)
    n=size(x);pow2=2
    if(n<2 .or. lambda<=0.0_dp)then;c=0.0_dp;return;end if
    a=gfact*0.920_dp*lambda*exp(-log(real(n,dp))/7.0_dp)/sqrt2
    b=gfact*0.912_dp*lambda*exp(-log(real(n,dp))/9.0_dp)/sqrt2
    ik=kernel_integral_moment(1,ker=wgaus)
    td=0.0_dp
    do i=1,n
      do j=i,n
        z=(x(i)-x(j))/b
        td=td+merge(1.0_dp,2.0_dp,i==j)*kernel_convolution6(z,wgaus)/(ik*ik)
      end do
    end do
    td=-td/real(n*(n-1),dp)
    c1=kernel_convolution4(0.0_dp,wgaus)
    c2=kernel_integral_moment(1,pow2,wgaus)
    c3=kernel_convolution(0.0_dp,wgaus)
    sa=compsda(x,a)
    if(td<=0.0_dp .or. sa<=0.0_dp)then;c=b;return;end if
    c=b*exp(log(c1*ik/(c2*c3*gfact**4)*sa/td)/7.0_dp)*sqrt2
  end function width_sj

  pure subroutine kde_criterion(x,h,criterion,ker,c,res)
    real(dp),intent(in)::x(:),h
    integer,intent(in)::criterion,ker
    real(dp),intent(in),optional::c
    real(dp),intent(out)::res(3)
    integer::i,j,n
    integer::pow2(1)
    real(dp)::degfree,dfd,pen,s,r0,r1,d0,d1,ik,wij,cc
    n=size(x);res=0.0_dp
    if(n<2 .or. h<=0.0_dp)then;res=huge(1.0_dp);return;end if
    ik=kernel_integral_moment(1,ker=ker);cc=0.0_dp;if(present(c))cc=c
    select case(criterion)
    case(kde_aic)
      pen=2.0_dp
      do i=1,n
        r0=0.0_dp;d0=0.0_dp
        do j=1,n
          s=(x(i)-x(j))/h
          r0=r0+kernel_weight(s,ker)
          d0=d0+s*s*kernel_derivative_over_u(s,ker)
        end do
        d0=-(d0+r0)/(real(n,dp)*h*h*ik)
        r0=r0/(real(n,dp)*h*ik)
        if(r0<=0.0_dp)then;res=huge(1.0_dp);return;end if
        res(1)=res(1)-2.0_dp*log(r0)+pen*kernel_weight(0.0_dp,ker)/(real(n,dp)*h*ik*r0)
        res(2)=res(2)-2.0_dp*d0/r0-pen*kernel_weight(0.0_dp,ker)/(real(n,dp)*h*ik*r0)*(d0/r0+1.0_dp/h)
      end do
    case(kde_ocv)
      do i=1,n
        r0=0.0_dp;d0=0.0_dp
        do j=1,n
          if(i==j)cycle
          s=(x(i)-x(j))/h
          r0=r0+kernel_weight(s,ker)
          d0=d0+s*s*kernel_derivative_over_u(s,ker)
        end do
        d0=-(d0+r0)/(real(n-1,dp)*h*h)
        r0=r0/(real(n-1,dp)*h)
        if(r0<=0.0_dp)then;res=huge(1.0_dp);return;end if
        res(1)=res(1)-log(r0);res(2)=res(2)-d0/r0
      end do
    case(kde_lscv)
      r0=0.0_dp;r1=0.0_dp;d0=0.0_dp;d1=0.0_dp;degfree=0.0_dp
      do i=1,n
        dfd=0.0_dp
        do j=1,n
          s=(x(i)-x(j))/h
          wij=kernel_weight(s,ker);dfd=dfd+wij
          r0=r0+kernel_convolution(s,ker)
          d0=d0-s*s*kernel_convolution1(s,ker)
          if(i/=j)then;r1=r1+wij;d1=d1-s*s*kernel_derivative_over_u(s,ker);end if
        end do
        if(dfd>0.0_dp)degfree=degfree+1.0_dp/dfd
      end do
      d1=d1-r1;d0=d0-r0
      res(1)=r0/(real(n*n,dp)*h*ik*ik)-2.0_dp*r1/(real(n*(n-1),dp)*h*ik)
      res(2)=d0/(real(n*n,dp)*h*h*ik*ik)-2.0_dp*d1/(real(n*(n-1),dp)*h*h*ik)
      res(3)=degfree
    case(kde_bcv)
      r0=0.0_dp;d0=0.0_dp
      do i=1,n
        do j=i+1,n
          s=(x(i)-x(j))/h
          r0=r0+2.0_dp*kernel_convolution4(s,ker)
          d0=d0+2.0_dp*s*kernel_convolution5(s,ker)
        end do
      end do
      d0=(-d0-r0)/(real(n*n,dp)*h*h*ik*ik);r0=r0/(real(n*n,dp)*h*ik*ik)
      pow2=2;d1=kernel_integral_moment(1,pow2,ker);r1=kernel_convolution(0.0_dp,ker)
      res(1)=(d1*d1*r0/4.0_dp+r1/(real(n,dp)*h))/(ik*ik)
      res(2)=(d1*d1*d0/4.0_dp-r1/(real(n,dp)*h*h))/(ik*ik)
    case(kde_sjpi,kde_gkk)
      if(criterion==kde_sjpi)then;s=cc*exp(5.0_dp*log(h)/7.0_dp)/sqrt2
      else;s=exp(log(real(n,dp))/10.0_dp)*h
      end if
      d0=compsda(x,s);res(1)=d0
      if(d0>0.0_dp)res(2)=exp(log(kernel_wikk(wgaus,0)/(d0*real(n,dp)))/5.0_dp)-h
    case default
      res=huge(1.0_dp)
    end select
  end subroutine kde_criterion

  pure real(dp) function solve_bandwidth(x,h0,h1,criterion,c,ker) result(besth)
    real(dp),intent(in)::x(:),h0,h1,c
    integer,intent(in)::criterion,ker
    real(dp)::h(7),d(7),r(7),res(3),minr,fact
    integer::i,nc,it,k
    k=10;minr=huge(1.0_dp);besth=h0;fact=1.00001_dp
    h=0.0_dp;d=0.0_dp;r=0.0_dp
    h(7)=h0;call kde_criterion(x,h(7),criterion,ker,c,res);r(7)=res(1);d(7)=res(2);nc=0
    do i=0,k-1
      h(6)=h(7);r(6)=r(7);d(6)=d(7)
      h(7)=h0*exp(real(i+1,dp)*log(h1/h0)/real(k,dp))
      call kde_criterion(x,h(7),criterion,ker,c,res);r(7)=res(1);d(7)=res(2)
      if(d(6)*d(7)<0.0_dp)then
        h(3)=h(6);d(3)=d(6);r(3)=r(6);h(1)=h(3);d(1)=d(3);r(1)=r(3)
        h(4)=h(7);d(4)=d(7);r(4)=r(7);h(2)=h(4);d(2)=d(4);r(2)=r(4)
        do it=1,100
          if(.not.(h(4)>fact*h(3) .or. h(3)>fact*h(4)))exit
          h(5)=h(4)-d(4)*(h(4)-h(3))/(d(4)-d(3))
          if(h(5)<h(1) .or. h(5)>h(2))h(5)=0.5_dp*(h(1)+h(2))
          call kde_criterion(x,h(5),criterion,ker,c,res);r(5)=res(1);d(5)=res(2)
          h(3)=h(4);d(3)=d(4);r(3)=r(4);h(4)=h(5);d(4)=d(5);r(4)=r(5)
          if(d(5)*d(1)>0.0_dp)then;h(1)=h(5);d(1)=d(5);r(1)=r(5)
          else;h(2)=h(5);d(2)=d(5);r(2)=r(5)
          end if
        end do
        if(criterion>=kde_bcv)then;besth=h(5);return;end if
        if(r(5)<=minr)then;minr=r(5);besth=h(5);end if
        nc=nc+1
      end if
    end do
    if(nc==0)besth=merge(h0,h1,r(6)<r(7))
  end function solve_bandwidth

  pure subroutine kde_select(x,h0,h1,methods,ker,bandwidths)
    real(dp),intent(in)::x(:),h0,h1
    integer,intent(in)::methods(:),ker
    real(dp),intent(out)::bandwidths(:)
    real(dp),allocatable::xs(:)
    real(dp)::scale,c
    integer::n,k,i
    n=size(x);allocate(xs(n));xs=x;call sort_values(xs)
    k=max(1,n/4)
    scale=xs(max(1,n+1-k))-xs(k)
    if(scale<=0.0_dp)scale=maxval(xs)-minval(xs)
    c=width_sj(x,scale)
    do i=1,size(methods)
      bandwidths(i)=solve_bandwidth(x,h0,h1,methods(i),c,ker)
    end do
  end subroutine kde_select

end module locfit_bandwidth
