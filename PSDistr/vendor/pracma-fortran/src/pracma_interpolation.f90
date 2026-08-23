! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_interpolation
   use pracma_kinds, only : dp
   use pracma_status, only : pracma_ok, pracma_invalid_argument, pracma_singular
   use pracma_types, only : pchip_result
   use pracma_linalg, only : solve_linear
   implicit none
   private
   public :: interp1, interp1_linear, interp1_nearest, interp1_spline, interp1_pchip
   public :: interp2, neville, newton_interp, lagrange_interp, barycentric, barycentric_weights
   public :: barylag, barylag2d, cubicspline, pchip, akima, ppval, mkpp, deval
   public :: ratinterp, lebesgue_function, spinterp, cutpoints, piecewise_linear

contains

   function interp1(x,y,xi,method,extrapolate) result(yi)
      real(dp),intent(in)::x(:),y(:),xi(:)
      character(len=*),intent(in),optional::method
      logical,intent(in),optional::extrapolate
      real(dp),allocatable::yi(:)
      character(len=16)::m
      m='linear'; if(present(method))m=adjustl(method)
      select case(trim(m))
      case('nearest'); yi=interp1_nearest(x,y,xi,extrapolate)
      case('spline','cubic'); yi=interp1_spline(x,y,xi,extrapolate)
      case('pchip'); yi=interp1_pchip(x,y,xi,extrapolate)
      case default; yi=interp1_linear(x,y,xi,extrapolate)
      end select
   end function interp1

   function interp1_linear(x,y,xi,extrapolate) result(yi)
      real(dp),intent(in)::x(:),y(:),xi(:)
      logical,intent(in),optional::extrapolate
      real(dp),allocatable::yi(:)
      logical::ext
      integer::i,j,n
      ext=.false.; if(present(extrapolate))ext=extrapolate
      n=size(x); allocate(yi(size(xi))); yi=0.0_dp
      if(n<2 .or. size(y)/=n)then; yi=huge(1.0_dp); return; end if
      do i=1,size(xi)
         if(xi(i)<x(1))then
            if(ext)then; j=1; else; yi(i)=huge(1.0_dp); cycle; end if
         else if(xi(i)>x(n))then
            if(ext)then; j=n-1; else; yi(i)=huge(1.0_dp); cycle; end if
         else
            j=locate_interval(x,xi(i))
         end if
         yi(i)=y(j)+(xi(i)-x(j))*(y(j+1)-y(j))/(x(j+1)-x(j))
      end do
   end function interp1_linear

   function interp1_nearest(x,y,xi,extrapolate) result(yi)
      real(dp),intent(in)::x(:),y(:),xi(:)
      logical,intent(in),optional::extrapolate
      real(dp),allocatable::yi(:)
      logical::ext
      integer::i,j,n
      ext=.false.; if(present(extrapolate))ext=extrapolate
      n=size(x); allocate(yi(size(xi))); yi=0.0_dp
      if(n<1 .or. size(y)/=n)then; yi=huge(1.0_dp); return; end if
      do i=1,size(xi)
         if((xi(i)<x(1).or.xi(i)>x(n)).and..not.ext)then; yi(i)=huge(1.0_dp); cycle; end if
         if(xi(i)<=x(1))then; j=1
         else if(xi(i)>=x(n))then; j=n
         else
            j=locate_interval(x,xi(i))
            if(abs(xi(i)-x(j+1))<abs(xi(i)-x(j)))j=j+1
         end if
         yi(i)=y(j)
      end do
   end function interp1_nearest

   function interp1_spline(x,y,xi,extrapolate) result(yi)
      real(dp),intent(in)::x(:),y(:),xi(:)
      logical,intent(in),optional::extrapolate
      real(dp),allocatable::yi(:)
      type(pchip_result)::pp
      pp=cubicspline(x,y)
      yi=ppval(pp,xi,extrapolate)
   end function interp1_spline

   function interp1_pchip(x,y,xi,extrapolate) result(yi)
      real(dp),intent(in)::x(:),y(:),xi(:)
      logical,intent(in),optional::extrapolate
      real(dp),allocatable::yi(:)
      type(pchip_result)::pp
      pp=pchip(x,y)
      yi=ppval(pp,xi,extrapolate)
   end function interp1_pchip

   function interp2(x,y,z,xi,yi,extrapolate) result(zi)
      real(dp),intent(in)::x(:),y(:),z(:,:),xi(:),yi(:)
      logical,intent(in),optional::extrapolate
      real(dp),allocatable::zi(:)
      logical::ext
      integer::k,i,j,nx,ny
      real(dp)::tx,ty
      ext=.false.; if(present(extrapolate))ext=extrapolate
      nx=size(x); ny=size(y); allocate(zi(size(xi))); zi=huge(1.0_dp)
      if(size(yi)/=size(xi).or.size(z,1)/=ny.or.size(z,2)/=nx.or.nx<2.or.ny<2)return
      do k=1,size(xi)
         if(.not.ext.and.(xi(k)<x(1).or.xi(k)>x(nx).or.yi(k)<y(1).or.yi(k)>y(ny)))cycle
         i=locate_clamped(x,xi(k)); j=locate_clamped(y,yi(k))
         tx=(xi(k)-x(i))/(x(i+1)-x(i)); ty=(yi(k)-y(j))/(y(j+1)-y(j))
         zi(k)=(1-tx)*(1-ty)*z(j,i)+tx*(1-ty)*z(j,i+1)+(1-tx)*ty*z(j+1,i)+tx*ty*z(j+1,i+1)
      end do
   end function interp2

   function neville(x,y,xi) result(v)
      real(dp),intent(in)::x(:),y(:),xi
      real(dp)::v
      real(dp),allocatable::q(:)
      integer::i,j,n
      n=size(x); allocate(q(n)); q=y
      do j=1,n-1
         do i=1,n-j
            q(i)=((xi-x(i+j))*q(i)+(x(i)-xi)*q(i+1))/(x(i)-x(i+j))
         end do
      end do
      v=q(1)
   end function neville

   function newton_interp(x,y,xi,coefficients) result(v)
      real(dp),intent(in)::x(:),y(:),xi
      real(dp),allocatable,intent(out),optional::coefficients(:)
      real(dp)::v
      real(dp),allocatable::c(:)
      integer::i,j,n
      n=size(x); allocate(c(n)); c=y
      do j=2,n
         do i=n,j,-1
            c(i)=(c(i)-c(i-1))/(x(i)-x(i-j+1))
         end do
      end do
      v=c(n)
      do i=n-1,1,-1; v=c(i)+(xi-x(i))*v; end do
      if(present(coefficients))coefficients=c
   end function newton_interp

   function lagrange_interp(x,y,xi) result(v)
      real(dp),intent(in)::x(:),y(:),xi
      real(dp)::v,l
      integer::i,j,n
      n=size(x); v=0.0_dp
      do i=1,n
         l=1.0_dp
         do j=1,n
            if(j/=i)l=l*(xi-x(j))/(x(i)-x(j))
         end do
         v=v+y(i)*l
      end do
   end function lagrange_interp

   function barycentric_weights(x) result(w)
      real(dp),intent(in)::x(:)
      real(dp),allocatable::w(:)
      integer::i,j,n
      n=size(x); allocate(w(n)); w=1.0_dp
      do i=1,n
         do j=1,n
            if(j/=i)w(i)=w(i)/(x(i)-x(j))
         end do
      end do
   end function barycentric_weights

   function barycentric(x,y,xi,w) result(v)
      real(dp),intent(in)::x(:),y(:),xi
      real(dp),intent(in),optional::w(:)
      real(dp)::v,num,den,t
      real(dp),allocatable::ww(:)
      integer::i
      if(present(w))then; ww=w; else; ww=barycentric_weights(x); end if
      do i=1,size(x)
         if(abs(xi-x(i))<=epsilon(1.0_dp)*max(1.0_dp,abs(x(i))))then; v=y(i); return; end if
      end do
      num=0.0_dp; den=0.0_dp
      do i=1,size(x); t=ww(i)/(xi-x(i)); num=num+t*y(i); den=den+t; end do
      v=num/den
   end function barycentric

   function barylag(x,y,xi) result(yi)
      real(dp),intent(in)::x(:),y(:),xi(:)
      real(dp),allocatable::yi(:),w(:)
      integer::i
      w=barycentric_weights(x); allocate(yi(size(xi)))
      do i=1,size(xi); yi(i)=barycentric(x,y,xi(i),w); end do
   end function barylag

   function barylag2d(x,y,z,xi,yi) result(zi)
      real(dp),intent(in)::x(:),y(:),z(:,:),xi(:),yi(:)
      real(dp),allocatable::zi(:),tmp(:),wx(:),wy(:)
      integer::k,j
      wx=barycentric_weights(x); wy=barycentric_weights(y)
      allocate(zi(size(xi)),tmp(size(y)))
      do k=1,size(xi)
         do j=1,size(y); tmp(j)=barycentric(x,z(j,:),xi(k),wx); end do
         zi(k)=barycentric(y,tmp,yi(k),wy)
      end do
   end function barylag2d

   function cubicspline(x,y) result(pp)
      real(dp),intent(in)::x(:),y(:)
      type(pchip_result)::pp
      integer::n,i
      real(dp),allocatable::h(:),a(:,:),rhs(:),m(:)
      n=size(x)
      if(n<2.or.size(y)/=n.or.any(x(2:)-x(:n-1)<=0.0_dp))then; pp%status=pracma_invalid_argument; return; end if
      allocate(pp%breaks(n),pp%coefficients(n-1,4),h(n-1)); pp%breaks=x; h=x(2:)-x(:n-1)
      allocate(a(n,n),rhs(n),m(n)); a=0.0_dp; rhs=0.0_dp; a(1,1)=1.0_dp; a(n,n)=1.0_dp
      do i=2,n-1
         a(i,i-1)=h(i-1); a(i,i)=2.0_dp*(h(i-1)+h(i)); a(i,i+1)=h(i)
         rhs(i)=6.0_dp*((y(i+1)-y(i))/h(i)-(y(i)-y(i-1))/h(i-1))
      end do
      call solve_linear(a,rhs,m)
      do i=1,n-1
         pp%coefficients(i,1)=(m(i+1)-m(i))/(6.0_dp*h(i))
         pp%coefficients(i,2)=0.5_dp*m(i)
         pp%coefficients(i,3)=(y(i+1)-y(i))/h(i)-h(i)*(2.0_dp*m(i)+m(i+1))/6.0_dp
         pp%coefficients(i,4)=y(i)
      end do
      pp%status=pracma_ok
   end function cubicspline

   function pchip(x,y) result(pp)
      real(dp),intent(in)::x(:),y(:)
      type(pchip_result)::pp
      integer::n,i
      real(dp),allocatable::h(:),delta(:),d(:)
      real(dp)::w1,w2
      n=size(x)
      if(n<2.or.size(y)/=n.or.any(x(2:)-x(:n-1)<=0.0_dp))then; pp%status=pracma_invalid_argument; return; end if
      allocate(h(n-1),delta(n-1),d(n),pp%breaks(n),pp%coefficients(n-1,4)); pp%breaks=x
      h=x(2:)-x(:n-1); delta=(y(2:)-y(:n-1))/h
      if(n==2)then
         d=delta(1)
      else
         d(1)=endpoint_slope(h(1),h(2),delta(1),delta(2))
         d(n)=endpoint_slope(h(n-1),h(n-2),delta(n-1),delta(n-2))
         do i=2,n-1
            if(delta(i-1)*delta(i)<=0.0_dp)then; d(i)=0.0_dp
            else
               w1=2.0_dp*h(i)+h(i-1); w2=h(i)+2.0_dp*h(i-1)
               d(i)=(w1+w2)/(w1/delta(i-1)+w2/delta(i))
            end if
         end do
      end if
      do i=1,n-1
         pp%coefficients(i,1)=(d(i)+d(i+1)-2.0_dp*delta(i))/(h(i)*h(i))
         pp%coefficients(i,2)=(3.0_dp*delta(i)-2.0_dp*d(i)-d(i+1))/h(i)
         pp%coefficients(i,3)=d(i); pp%coefficients(i,4)=y(i)
      end do
      pp%status=pracma_ok
   end function pchip

   function akima(x,y) result(pp)
      real(dp),intent(in)::x(:),y(:)
      type(pchip_result)::pp
      integer::n,i
      real(dp),allocatable::m(:),d(:),h(:)
      real(dp)::w1,w2
      n=size(x)
      if(n<5.or.size(y)/=n.or.any(x(2:)-x(:n-1)<=0.0_dp))then
         pp=pchip(x,y); return
      end if
      allocate(m(n+3),d(n),h(n-1),pp%breaks(n),pp%coefficients(n-1,4)); h=x(2:)-x(:n-1)
      m(3:n+1)=(y(2:)-y(:n-1))/h
      m(2)=2*m(3)-m(4); m(1)=2*m(2)-m(3); m(n+2)=2*m(n+1)-m(n); m(n+3)=2*m(n+2)-m(n+1)
      do i=1,n
         w1=abs(m(i+3)-m(i+2)); w2=abs(m(i+1)-m(i))
         if(w1+w2>tiny(1.0_dp))then; d(i)=(w1*m(i+1)+w2*m(i+2))/(w1+w2)
         else; d(i)=0.5_dp*(m(i+1)+m(i+2)); end if
      end do
      pp%breaks=x
      do i=1,n-1
         pp%coefficients(i,1)=(d(i)+d(i+1)-2.0_dp*m(i+2))/(h(i)*h(i))
         pp%coefficients(i,2)=(3.0_dp*m(i+2)-2*d(i)-d(i+1))/h(i)
         pp%coefficients(i,3)=d(i); pp%coefficients(i,4)=y(i)
      end do
      pp%status=pracma_ok
   end function akima

   function mkpp(breaks,coefficients) result(pp)
      real(dp),intent(in)::breaks(:),coefficients(:,:)
      type(pchip_result)::pp
      if(size(coefficients,1)/=size(breaks)-1)then; pp%status=pracma_invalid_argument; return; end if
      pp%breaks=breaks; pp%coefficients=coefficients; pp%status=pracma_ok
   end function mkpp

   function ppval(pp,x,extrapolate) result(y)
      type(pchip_result),intent(in)::pp
      real(dp),intent(in)::x(:)
      logical,intent(in),optional::extrapolate
      real(dp),allocatable::y(:)
      logical::ext
      integer::i,j,k,nc,n
      real(dp)::dx,v
      ext=.false.; if(present(extrapolate))ext=extrapolate
      allocate(y(size(x))); y=huge(1.0_dp)
      if(.not.allocated(pp%breaks).or..not.allocated(pp%coefficients))return
      n=size(pp%breaks); nc=size(pp%coefficients,2)
      do i=1,size(x)
         if(.not.ext.and.(x(i)<pp%breaks(1).or.x(i)>pp%breaks(n)))cycle
         j=locate_clamped(pp%breaks,x(i)); dx=x(i)-pp%breaks(j); v=pp%coefficients(j,1)
         do k=2,nc; v=v*dx+pp%coefficients(j,k); end do
         y(i)=v
      end do
   end function ppval

   function deval(pp,x) result(y)
      type(pchip_result),intent(in)::pp
      real(dp),intent(in)::x(:)
      real(dp),allocatable::y(:)
      y=ppval(pp,x,.true.)
   end function deval

   function ratinterp(x,y,xi) result(v)
      real(dp),intent(in)::x(:),y(:),xi
      real(dp)::v
      real(dp),allocatable::c(:),d(:)
      real(dp)::w,t,dd
      integer::n,i,m,ns
      n=size(x); allocate(c(n),d(n)); ns=1
      dd=abs(xi-x(1))
      do i=1,n
         t=abs(xi-x(i)); if(t==0.0_dp)then; v=y(i); return; end if
         if(t<dd)then; ns=i; dd=t; end if
         c(i)=y(i); d(i)=y(i)+tiny(1.0_dp)
      end do
      v=y(ns); ns=ns-1
      do m=1,n-1
         do i=1,n-m
            w=c(i+1)-d(i); t=(x(i)-xi)*d(i)/(x(i+m)-xi); dd=t-c(i+1)
            if(abs(dd)<=tiny(1.0_dp))then; v=huge(1.0_dp); return; end if
            dd=w/dd; d(i)=c(i+1)*dd; c(i)=t*dd
         end do
         if(2*ns<n-m)then; v=v+c(ns+1)
         else; v=v+d(ns); ns=ns-1; end if
      end do
   end function ratinterp

   function lebesgue_function(x,xi) result(l)
      real(dp),intent(in)::x(:),xi(:)
      real(dp),allocatable::l(:),w(:)
      real(dp)::den,t
      integer::i,j
      w=barycentric_weights(x); allocate(l(size(xi)))
      do i=1,size(xi)
         if(any(abs(xi(i)-x)<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(x)))))then; l(i)=1.0_dp; cycle; end if
         den=sum(w/(xi(i)-x)); l(i)=0.0_dp
         do j=1,size(x); t=w(j)/(xi(i)-x(j))/den; l(i)=l(i)+abs(t); end do
      end do
   end function lebesgue_function

   function spinterp(x,y,xi,method) result(yi)
      real(dp),intent(in)::x(:),y(:),xi(:)
      character(len=*),intent(in),optional::method
      real(dp),allocatable::yi(:)
      yi=interp1(x,y,xi,method,.true.)
   end function spinterp

   function cutpoints(x,nparts) result(cuts)
      real(dp),intent(in)::x(:)
      integer,intent(in)::nparts
      real(dp),allocatable::cuts(:),s(:)
      integer::i,n,k
      n=size(x); if(nparts<1.or.n<1)then; allocate(cuts(0)); return; end if
      s=x; call sort_inplace(s); allocate(cuts(max(0,nparts-1)))
      do i=1,nparts-1; k=max(1,min(n,int(real(i*n,dp)/real(nparts,dp)))); cuts(i)=s(k); end do
   end function cutpoints

   function piecewise_linear(x,breaks,values) result(y)
      real(dp),intent(in)::x(:),breaks(:),values(:)
      real(dp),allocatable::y(:)
      y=interp1_linear(breaks,values,x,.true.)
   end function piecewise_linear

   pure integer function locate_interval(x,v) result(j)
      real(dp),intent(in)::x(:),v
      integer::lo,hi,mid,n
      n=size(x); if(v<=x(1))then; j=1; return; end if
      if(v>=x(n))then; j=n-1; return; end if
      lo=1; hi=n
      do while(hi-lo>1); mid=(lo+hi)/2; if(v>=x(mid))then; lo=mid; else; hi=mid; end if; end do
      j=lo
   end function locate_interval

   pure integer function locate_clamped(x,v) result(j)
      real(dp),intent(in)::x(:),v
      j=locate_interval(x,v)
   end function locate_clamped

   pure real(dp) function endpoint_slope(h1,h2,d1,d2) result(d)
      real(dp),intent(in)::h1,h2,d1,d2
      d=((2*h1+h2)*d1-h1*d2)/(h1+h2)
      if(d*d1<=0.0_dp)d=0.0_dp
      if(d1*d2<0.0_dp.and.abs(d)>3.0_dp*abs(d1))d=3.0_dp*d1
   end function endpoint_slope

   subroutine sort_inplace(a)
      real(dp),intent(inout)::a(:)
      integer::i,j
      real(dp)::t
      do i=2,size(a); t=a(i); j=i-1; do while(j>=1); if(a(j)<=t)exit; a(j+1)=a(j); j=j-1; end do; a(j+1)=t; end do
   end subroutine sort_inplace

end module pracma_interpolation
