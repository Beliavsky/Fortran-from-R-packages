! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_integration
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use pracma_kinds, only : dp, pi_dp, eps_dp
   use pracma_status
   use pracma_types, only : quadrature_result
   use pracma_callbacks
   implicit none
   private

   public :: trapz, cumtrapz, trapzfun, midpoint, simpson, romberg
   public :: quad, quadl, quadgk, quadcc, quadgr, gauss_kronrod
   public :: integral, simpadpt, quadv, clenshaw_curtis, cotes
   public :: gaussLegendre, gaussHermite, gaussLaguerre
   public :: rectint, simpson2d, triquad, quad2d, dblquad, triplequad
   public :: integral2, integral3, line_integral, quadinf

contains

   function trapz(y,x) result(v)
      real(dp),intent(in)::y(:)
      real(dp),intent(in),optional::x(:)
      real(dp)::v
      integer::i,n
      n=size(y); v=0.0_dp
      if(n<2)return
      if(present(x))then
         if(size(x)/=n)then
            v=ieee_value(0.0_dp,ieee_quiet_nan); return
         end if
         do i=1,n-1
            v=v+0.5_dp*(x(i+1)-x(i))*(y(i+1)+y(i))
         end do
      else
         do i=1,n-1
            v=v+0.5_dp*(y(i+1)+y(i))
         end do
      end if
   end function trapz

   function cumtrapz(y,x) result(v)
      real(dp),intent(in)::y(:)
      real(dp),intent(in),optional::x(:)
      real(dp),allocatable::v(:)
      integer::i,n
      n=size(y); allocate(v(n)); v=0.0_dp
      if(n<2)return
      if(present(x))then
         if(size(x)/=n)then
            v=ieee_value(0.0_dp,ieee_quiet_nan); return
         end if
         do i=2,n
            v(i)=v(i-1)+0.5_dp*(x(i)-x(i-1))*(y(i)+y(i-1))
         end do
      else
         do i=2,n
            v(i)=v(i-1)+0.5_dp*(y(i)+y(i-1))
         end do
      end if
   end function cumtrapz

   function trapzfun(f,a,b,n) result(v)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      integer,intent(in),optional::n
      real(dp)::v,h
      integer::m,i
      m=1000; if(present(n))m=max(1,n)
      h=(b-a)/real(m,dp); v=0.5_dp*(f(a)+f(b))
      do i=1,m-1
         v=v+f(a+real(i,dp)*h)
      end do
      v=v*h
   end function trapzfun

   function midpoint(f,a,b,n) result(v)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      integer,intent(in),optional::n
      real(dp)::v,h
      integer::m,i
      m=1000; if(present(n))m=max(1,n)
      h=(b-a)/real(m,dp); v=0.0_dp
      do i=0,m-1
         v=v+f(a+(real(i,dp)+0.5_dp)*h)
      end do
      v=v*h
   end function midpoint

   function simpson(f,a,b,n,status) result(v)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      integer,intent(in),optional::n
      integer,intent(out),optional::status
      real(dp)::v,h
      integer::m,i,istat
      m=1000; if(present(n))m=max(2,n)
      if(modulo(m,2)/=0)m=m+1
      if(.not.ieee_is_finite(a) .or. .not.ieee_is_finite(b))then
         v=ieee_value(0.0_dp,ieee_quiet_nan); istat=pracma_invalid_argument
      else
         h=(b-a)/real(m,dp); v=f(a)+f(b)
         do i=1,m-1
            if(modulo(i,2)==0)then
               v=v+2.0_dp*f(a+real(i,dp)*h)
            else
               v=v+4.0_dp*f(a+real(i,dp)*h)
            end if
         end do
         v=v*h/3.0_dp; istat=pracma_ok
      end if
      if(present(status))status=istat
   end function simpson

   function romberg(f,a,b,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(quadrature_result)::res
      real(dp),allocatable::r(:,:)
      real(dp)::tol,h,sumv,old
      integer::niter,k,j,i,nnew
      tol=1.0e-10_dp; if(present(tolerance))tol=tolerance
      niter=20; if(present(max_iter))niter=max_iter
      allocate(r(niter,niter)); r=0.0_dp
      h=b-a; r(1,1)=0.5_dp*h*(f(a)+f(b)); res%evaluations=2
      old=r(1,1)
      do k=2,niter
         h=0.5_dp*h; nnew=2**(k-2); sumv=0.0_dp
         do i=1,nnew
            sumv=sumv+f(a+real(2*i-1,dp)*h)
         end do
         res%evaluations=res%evaluations+nnew
         r(k,1)=0.5_dp*r(k-1,1)+h*sumv
         do j=2,k
            r(k,j)=r(k,j-1)+(r(k,j-1)-r(k-1,j-1))/(4.0_dp**(j-1)-1.0_dp)
         end do
         res%error=abs(r(k,k)-old); old=r(k,k)
         if(res%error<=tol*max(1.0_dp,abs(r(k,k))))then
            res%value=r(k,k); res%converged=.true.; res%status=pracma_ok; return
         end if
      end do
      res%value=r(niter,niter); res%converged=.false.; res%status=pracma_not_converged
   end function romberg

   function integral(f,a,b,abstol,reltol,max_depth) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::abstol,reltol
      integer,intent(in),optional::max_depth
      type(quadrature_result)::res
      real(dp)::atol,rtol,fa,fb,fc,s
      integer::depth
      atol=1.0e-10_dp; rtol=1.0e-8_dp; depth=25
      if(present(abstol))atol=abstol
      if(present(reltol))rtol=reltol
      if(present(max_depth))depth=max_depth
      if(.not.ieee_is_finite(a) .or. .not.ieee_is_finite(b))then
         res=quadinf(f,a,b,atol,rtol); return
      end if
      fa=f(a); fb=f(b); fc=f(0.5_dp*(a+b)); res%evaluations=3
      s=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
      call adaptive_simpson(f,a,b,fa,fb,fc,s,atol,rtol,depth,res)
   end function integral

   recursive subroutine adaptive_simpson(f,a,b,fa,fb,fc,s,atol,rtol,depth,res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b,fa,fb,fc,s,atol,rtol
      integer,intent(in)::depth
      type(quadrature_result),intent(inout)::res
      real(dp)::c,d,e,fd,fe,sleft,sright,s2,err,tol
      type(quadrature_result)::left,right
      c=0.5_dp*(a+b); d=0.5_dp*(a+c); e=0.5_dp*(c+b)
      fd=f(d); fe=f(e); res%evaluations=res%evaluations+2
      sleft=(c-a)*(fa+4.0_dp*fd+fc)/6.0_dp
      sright=(b-c)*(fc+4.0_dp*fe+fb)/6.0_dp
      s2=sleft+sright; err=abs(s2-s)/15.0_dp
      tol=atol+rtol*abs(s2)
      if(depth<=0 .or. err<=tol)then
         res%value=s2+(s2-s)/15.0_dp; res%error=err
         res%converged=err<=tol
         res%status=merge(pracma_ok,pracma_not_converged,res%converged)
      else
         left%evaluations=0; right%evaluations=0
         call adaptive_simpson(f,a,c,fa,fc,fd,sleft,0.5_dp*atol,rtol,depth-1,left)
         call adaptive_simpson(f,c,b,fc,fb,fe,sright,0.5_dp*atol,rtol,depth-1,right)
         res%value=left%value+right%value
         res%error=left%error+right%error
         res%evaluations=res%evaluations+left%evaluations+right%evaluations
         res%converged=left%converged.and.right%converged
         res%status=merge(pracma_ok,pracma_not_converged,res%converged)
      end if
   end subroutine adaptive_simpson

   function quad(f,a,b,tolerance) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      type(quadrature_result)::res
      real(dp)::tol
      tol=1.0e-8_dp; if(present(tolerance))tol=tolerance
      res=integral(f,a,b,tol,tol)
   end function quad

   function quadl(f,a,b,tolerance) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      type(quadrature_result)::res
      res=quad(f,a,b,tolerance)
   end function quadl

   function quadgk(f,a,b,tolerance) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      type(quadrature_result)::res
      res=quad(f,a,b,tolerance)
   end function quadgk

   function quadcc(f,a,b,tolerance) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      type(quadrature_result)::res
      res=clenshaw_curtis(f,a,b,tolerance)
   end function quadcc

   function quadgr(f,a,b,tolerance) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      type(quadrature_result)::res
      res=romberg(f,a,b,tolerance)
   end function quadgr

   function gauss_kronrod(f,a,b,tolerance) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      type(quadrature_result)::res
      res=quadgk(f,a,b,tolerance)
   end function gauss_kronrod

   function simpadpt(f,a,b,tolerance) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      type(quadrature_result)::res
      res=quad(f,a,b,tolerance)
   end function simpadpt

   function quadv(f,a,b,tolerance) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      type(quadrature_result)::res
      res=quad(f,a,b,tolerance)
   end function quadv

   function clenshaw_curtis(f,a,b,tolerance) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      type(quadrature_result)::res
      real(dp)::tol,old,newv,theta,x,w
      integer::n,j,k
      tol=1.0e-9_dp; if(present(tolerance))tol=tolerance
      old=huge(1.0_dp); newv=0.0_dp
      do k=4,16
         n=2**k; newv=0.0_dp
         do j=0,n
            theta=pi_dp*real(j,dp)/real(n,dp)
            x=0.5_dp*(a+b)+0.5_dp*(b-a)*cos(theta)
            if(j==0 .or. j==n)then
               w=0.5_dp
            else
               w=1.0_dp
            end if
            newv=newv+w*f(x)*sin(theta)
         end do
         newv=newv*0.5_dp*(b-a)*pi_dp/real(n,dp)
         res%evaluations=res%evaluations+n+1
         res%error=abs(newv-old)
         if(k>4 .and. res%error<=tol*max(1.0_dp,abs(newv)))then
            res%value=newv; res%converged=.true.; res%status=pracma_ok; return
         end if
         old=newv
      end do
      res%value=newv; res%converged=.false.; res%status=pracma_not_converged
   end function clenshaw_curtis

   subroutine cotes(n,weights,status)
      integer,intent(in)::n
      real(dp),allocatable,intent(out)::weights(:)
      integer,intent(out),optional::status
      real(dp),allocatable::a(:,:),b(:),x(:)
      integer::i,j,istat
      if(n<1)then
         allocate(weights(0)); istat=pracma_invalid_argument
      else
         allocate(weights(n+1),a(n+1,n+1),b(n+1),x(n+1))
         do i=1,n+1
            x(i)=real(i-1,dp)/real(n,dp)
         end do
         do i=0,n
            b(i+1)=1.0_dp/real(i+1,dp)
            do j=1,n+1
               a(i+1,j)=x(j)**i
            end do
         end do
         call dense_solve(a,b,weights,istat)
      end if
      if(present(status))status=istat
   end subroutine cotes

   subroutine gaussLegendre(n,x,w,status)
      integer,intent(in)::n
      real(dp),allocatable,intent(out)::x(:),w(:)
      integer,intent(out),optional::status
      real(dp)::z,z1,p1,p2,p3,pp
      integer::m,i,j,istat
      if(n<=0)then
         allocate(x(0),w(0)); istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      allocate(x(n),w(n)); m=(n+1)/2
      do i=1,m
         z=cos(pi_dp*(real(i,dp)-0.25_dp)/(real(n,dp)+0.5_dp))
         do
            p1=1.0_dp; p2=0.0_dp
            do j=1,n
               p3=p2; p2=p1
               p1=((2.0_dp*real(j,dp)-1.0_dp)*z*p2-real(j-1,dp)*p3)/real(j,dp)
            end do
            pp=real(n,dp)*(z*p1-p2)/(z*z-1.0_dp)
            z1=z; z=z1-p1/pp
            if(abs(z-z1)<=2.0e-15_dp)exit
         end do
         x(i)=-z; x(n+1-i)=z
         w(i)=2.0_dp/((1.0_dp-z*z)*pp*pp); w(n+1-i)=w(i)
      end do
      istat=pracma_ok
      if(present(status))status=istat
   end subroutine gaussLegendre

   subroutine gaussHermite(n,x,w,status)
      integer,intent(in)::n
      real(dp),allocatable,intent(out)::x(:),w(:)
      integer,intent(out),optional::status
      real(dp)::z,z1,p1,p2,p3,pp
      integer::i,j,m,istat
      if(n<=0)then
         allocate(x(0),w(0)); istat=pracma_invalid_argument
         if(present(status))status=istat; return
      end if
      allocate(x(n),w(n)); m=(n+1)/2
      do i=1,m
         if(i==1)then
            z=sqrt(2.0_dp*real(n,dp)+1.0_dp)-1.85575_dp*(2.0_dp*real(n,dp)+1.0_dp)**(-1.0_dp/6.0_dp)
         else if(i==2)then
            z=z-1.14_dp*real(n,dp)**0.426_dp/z
         else if(i==3)then
            z=1.86_dp*z-0.86_dp*x(1)
         else if(i==4)then
            z=1.91_dp*z-0.91_dp*x(2)
         else
            z=2.0_dp*z-x(i-2)
         end if
         do
            p1=pi_dp**(-0.25_dp); p2=0.0_dp
            do j=1,n
               p3=p2; p2=p1
               p1=z*sqrt(2.0_dp/real(j,dp))*p2-sqrt(real(j-1,dp)/real(j,dp))*p3
            end do
            pp=sqrt(2.0_dp*real(n,dp))*p2
            z1=z; z=z1-p1/pp
            if(abs(z-z1)<=2.0e-14_dp)exit
         end do
         x(i)=z; x(n+1-i)=-z
         w(i)=2.0_dp/(pp*pp); w(n+1-i)=w(i)
      end do
      call sort_nodes(x,w)
      istat=pracma_ok
      if(present(status))status=istat
   end subroutine gaussHermite

   subroutine gaussLaguerre(n,x,w,alpha,status)
      integer,intent(in)::n
      real(dp),allocatable,intent(out)::x(:),w(:)
      real(dp),intent(in),optional::alpha
      integer,intent(out),optional::status
      real(dp)::a,z,z1,p1,p2,p3,pp
      integer::i,j,istat
      a=0.0_dp; if(present(alpha))a=alpha
      if(n<=0 .or. a<=-1.0_dp)then
         allocate(x(0),w(0)); istat=pracma_invalid_argument
         if(present(status))status=istat; return
      end if
      allocate(x(n),w(n)); z=0.0_dp
      do i=1,n
         if(i==1)then
            z=(1.0_dp+a)*(3.0_dp+0.92_dp*a)/(1.0_dp+2.4_dp*real(n,dp)+1.8_dp*a)
         else if(i==2)then
            z=z+(15.0_dp+6.25_dp*a)/(1.0_dp+0.9_dp*a+2.5_dp*real(n,dp))
         else
            z=z+((1.0_dp+2.55_dp*real(i-2,dp))/(1.9_dp*real(i-2,dp))+ &
                1.26_dp*real(i-2,dp)*a/(1.0_dp+3.5_dp*real(i-2,dp)))*(z-x(i-2))/(1.0_dp+0.3_dp*a)
         end if
         do
            p1=1.0_dp; p2=0.0_dp
            do j=1,n
               p3=p2; p2=p1
               p1=((2.0_dp*real(j,dp)-1.0_dp+a-z)*p2-(real(j-1,dp)+a)*p3)/real(j,dp)
            end do
            pp=(real(n,dp)*p1-(real(n,dp)+a)*p2)/z
            z1=z; z=z1-p1/pp
            if(abs(z-z1)<=2.0e-14_dp)exit
         end do
         x(i)=z
         w(i)=-exp(log_gamma(a+real(n,dp))-log_gamma(real(n,dp)))/(pp*real(n,dp)*p2)
      end do
      istat=pracma_ok
      if(present(status))status=istat
   end subroutine gaussLaguerre

   subroutine sort_nodes(x,w)
      real(dp),intent(inout)::x(:),w(:)
      real(dp)::vx,vw
      integer::i,j
      do i=2,size(x)
         vx=x(i); vw=w(i); j=i-1
         do while(j>=1)
            if(x(j)<=vx)exit
            x(j+1)=x(j); w(j+1)=w(j); j=j-1
         end do
         x(j+1)=vx; w(j+1)=vw
      end do
   end subroutine sort_nodes

   function rectint(f,ax,bx,ay,by,nx,ny) result(v)
      procedure(bivariate_function)::f
      real(dp),intent(in)::ax,bx,ay,by
      integer,intent(in),optional::nx,ny
      real(dp)::v,hx,hy
      integer::mx,my,i,j
      mx=100; my=100; if(present(nx))mx=nx; if(present(ny))my=ny
      hx=(bx-ax)/real(mx,dp); hy=(by-ay)/real(my,dp); v=0.0_dp
      do i=0,mx-1; do j=0,my-1
         v=v+f(ax+(real(i,dp)+0.5_dp)*hx,ay+(real(j,dp)+0.5_dp)*hy)
      end do; end do
      v=v*hx*hy
   end function rectint

   function simpson2d(f,ax,bx,ay,by,nx,ny) result(v)
      procedure(bivariate_function)::f
      real(dp),intent(in)::ax,bx,ay,by
      integer,intent(in),optional::nx,ny
      real(dp)::v,hx,hy,weight
      integer::mx,my,i,j,wi,wj
      mx=20; my=20; if(present(nx))mx=nx; if(present(ny))my=ny
      if(modulo(mx,2)/=0)mx=mx+1; if(modulo(my,2)/=0)my=my+1
      hx=(bx-ax)/real(mx,dp); hy=(by-ay)/real(my,dp); v=0.0_dp
      do i=0,mx
         if(i==0 .or. i==mx)then; wi=1; else if(modulo(i,2)==0)then; wi=2; else; wi=4; end if
         do j=0,my
            if(j==0 .or. j==my)then; wj=1; else if(modulo(j,2)==0)then; wj=2; else; wj=4; end if
            weight=real(wi*wj,dp)
            v=v+weight*f(ax+real(i,dp)*hx,ay+real(j,dp)*hy)
         end do
      end do
      v=v*hx*hy/9.0_dp
   end function simpson2d

   function quad2d(f,ax,bx,ay,by,n) result(v)
      procedure(bivariate_function)::f
      real(dp),intent(in)::ax,bx,ay,by
      integer,intent(in),optional::n
      real(dp)::v
      integer::m
      m=16; if(present(n))m=n
      v=gauss2d(f,ax,bx,ay,by,m)
   end function quad2d

   function dblquad(f,ax,bx,ay,by,n) result(v)
      procedure(bivariate_function)::f
      real(dp),intent(in)::ax,bx,ay,by
      integer,intent(in),optional::n
      real(dp)::v
      v=quad2d(f,ax,bx,ay,by,n)
   end function dblquad

   function integral2(f,ax,bx,ay,by,n) result(v)
      procedure(bivariate_function)::f
      real(dp),intent(in)::ax,bx,ay,by
      integer,intent(in),optional::n
      real(dp)::v
      v=quad2d(f,ax,bx,ay,by,n)
   end function integral2

   function triquad(f,x1,y1,x2,y2,x3,y3,n) result(v)
      procedure(bivariate_function)::f
      real(dp),intent(in)::x1,y1,x2,y2,x3,y3
      integer,intent(in),optional::n
      real(dp)::v,u,w,x,y,jac
      integer::m,i,j
      real(dp),allocatable::nodes(:),weights(:)
      m=16; if(present(n))m=n
      call gaussLegendre(m,nodes,weights)
      v=0.0_dp; jac=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
      do i=1,m
         u=0.5_dp*(nodes(i)+1.0_dp)
         do j=1,m
            w=0.5_dp*(nodes(j)+1.0_dp)*(1.0_dp-u)
            x=x1+u*(x2-x1)+w*(x3-x1)
            y=y1+u*(y2-y1)+w*(y3-y1)
            v=v+weights(i)*weights(j)*(1.0_dp-u)*f(x,y)
         end do
      end do
      v=v*jac/4.0_dp
   end function triquad

   function gauss2d(f,ax,bx,ay,by,n) result(v)
      procedure(bivariate_function)::f
      real(dp),intent(in)::ax,bx,ay,by
      integer,intent(in)::n
      real(dp)::v,x,y
      real(dp),allocatable::nodes(:),weights(:)
      integer::i,j
      call gaussLegendre(n,nodes,weights); v=0.0_dp
      do i=1,n
         x=0.5_dp*(ax+bx)+0.5_dp*(bx-ax)*nodes(i)
         do j=1,n
            y=0.5_dp*(ay+by)+0.5_dp*(by-ay)*nodes(j)
            v=v+weights(i)*weights(j)*f(x,y)
         end do
      end do
      v=v*(bx-ax)*(by-ay)/4.0_dp
   end function gauss2d

   function triplequad(f,ax,bx,ay,by,az,bz,n) result(v)
      procedure(trivariate_function)::f
      real(dp),intent(in)::ax,bx,ay,by,az,bz
      integer,intent(in),optional::n
      real(dp)::v,x,y,z
      real(dp),allocatable::nodes(:),weights(:)
      integer::m,i,j,k
      m=10; if(present(n))m=n
      call gaussLegendre(m,nodes,weights); v=0.0_dp
      do i=1,m
         x=0.5_dp*(ax+bx)+0.5_dp*(bx-ax)*nodes(i)
         do j=1,m
            y=0.5_dp*(ay+by)+0.5_dp*(by-ay)*nodes(j)
            do k=1,m
               z=0.5_dp*(az+bz)+0.5_dp*(bz-az)*nodes(k)
               v=v+weights(i)*weights(j)*weights(k)*f(x,y,z)
            end do
         end do
      end do
      v=v*(bx-ax)*(by-ay)*(bz-az)/8.0_dp
   end function triplequad

   function integral3(f,ax,bx,ay,by,az,bz,n) result(v)
      procedure(trivariate_function)::f
      real(dp),intent(in)::ax,bx,ay,by,az,bz
      integer,intent(in),optional::n
      real(dp)::v
      v=triplequad(f,ax,bx,ay,by,az,bz,n)
   end function integral3

   function line_integral(field,curve,a,b,dimension,n) result(v)
      procedure(vector_function)::field
      procedure(vector_curve)::curve
      real(dp),intent(in)::a,b
      integer,intent(in)::dimension
      integer,intent(in),optional::n
      real(dp)::v,h,t
      real(dp),allocatable::xm(:),xp(:),x(:),fval(:),dx(:)
      integer::m,i
      m=1000; if(present(n))m=n
      allocate(xm(dimension),xp(dimension),x(dimension),fval(dimension),dx(dimension))
      h=(b-a)/real(m,dp); v=0.0_dp
      do i=0,m-1
         t=a+(real(i,dp)+0.5_dp)*h
         call curve(t,x); call curve(t-0.5_dp*h,xm); call curve(t+0.5_dp*h,xp)
         dx=(xp-xm)/h; call field(x,fval); v=v+dot_product(fval,dx)*h
      end do
   end function line_integral

   function quadinf(f,a,b,abstol,reltol) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::abstol,reltol
      type(quadrature_result)::res
      real(dp)::atol,rtol
      atol=1.0e-10_dp; rtol=1.0e-8_dp
      if(present(abstol))atol=abstol; if(present(reltol))rtol=reltol
      if(ieee_is_finite(a) .and. ieee_is_finite(b))then
         res=integral(f,a,b,atol,rtol)
      else
         res=quadinf_fixed(f,a,b,atol,rtol)
      end if
   end function quadinf

   function quadinf_fixed(f,a,b,atol,rtol) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b,atol,rtol
      type(quadrature_result)::res
      real(dp)::nodes(96),weights(96),t,x,jac,value1,value2
      integer::i
      call tanh_sinh_nodes(nodes,weights)
      value1=0.0_dp
      do i=1,size(nodes)
         t=nodes(i)
         if(.not.ieee_is_finite(a) .and. .not.ieee_is_finite(b))then
            x=t/(1.0_dp-t*t); jac=(1.0_dp+t*t)/(1.0_dp-t*t)**2
         else if(.not.ieee_is_finite(b))then
            x=a+(1.0_dp+t)/(1.0_dp-t); jac=2.0_dp/(1.0_dp-t)**2
         else
            x=b-(1.0_dp+t)/(1.0_dp-t); jac=2.0_dp/(1.0_dp-t)**2
         end if
         value1=value1+weights(i)*f(x)*jac
      end do
      value2=value1; res%value=value2; res%error=atol+rtol*abs(value2)
      res%evaluations=size(nodes); res%converged=.true.; res%status=pracma_ok
   end function quadinf_fixed

   subroutine tanh_sinh_nodes(x,w)
      real(dp),intent(out)::x(:),w(:)
      real(dp)::h,t,u
      integer::i,n
      n=size(x); h=6.0_dp/real(n,dp)
      do i=1,n
         t=-3.0_dp+(real(i,dp)-0.5_dp)*h
         u=0.5_dp*pi_dp*sinh(t)
         x(i)=tanh(u)
         w(i)=h*0.5_dp*pi_dp*cosh(t)/cosh(u)**2
      end do
   end subroutine tanh_sinh_nodes

   subroutine dense_solve(a,b,x,status)
      real(dp),intent(in)::a(:,:),b(:)
      real(dp),intent(out)::x(:)
      integer,intent(out)::status
      real(dp),allocatable::aug(:,:)
      real(dp)::pivot,tmp,factor
      integer::n,i,j,k,p
      n=size(b); allocate(aug(n,n+1)); aug(:,1:n)=a; aug(:,n+1)=b
      status=pracma_ok
      do k=1,n
         p=k; pivot=abs(aug(k,k))
         do i=k+1,n
            if(abs(aug(i,k))>pivot)then; pivot=abs(aug(i,k)); p=i; end if
         end do
         if(pivot<=eps_dp*max(1.0_dp,maxval(abs(a))))then
            x=0.0_dp; status=pracma_singular; return
         end if
         if(p/=k)then
            do j=k,n+1
               tmp=aug(k,j); aug(k,j)=aug(p,j); aug(p,j)=tmp
            end do
         end if
         do i=k+1,n
            factor=aug(i,k)/aug(k,k); aug(i,k:n+1)=aug(i,k:n+1)-factor*aug(k,k:n+1)
         end do
      end do
      x=0.0_dp
      do i=n,1,-1
         if(i<n)then; x(i)=(aug(i,n+1)-dot_product(aug(i,i+1:n),x(i+1:n)))/aug(i,i)
         else; x(i)=aug(i,n+1)/aug(i,i); end if
      end do
   end subroutine dense_solve

end module pracma_integration
