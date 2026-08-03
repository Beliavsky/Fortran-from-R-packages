! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_polynomial
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use pracma_kinds, only : dp, pi_dp, eps_dp
   use pracma_status
   use pracma_types, only : polynomial_division_result, polynomial_fit_result
   use pracma_linalg, only : mldivide, charpoly
   use pracma_basic, only : linspace
   implicit none
   private

   public :: polyval, polyvalm, polyadd, polymul, polypow, polydiv
   public :: polyder, polyint, polytrans, poly2str, roots, rootsmult
   public :: polyroots, compan, horner, hornerdefl, polyfit, polyfix
   public :: chebPoly, chebCoeff, chebApprox, legendre, laguerre
   public :: bernstein, bernsteinb, pade, Poly, polygcf, rationalfit
   public :: trigPoly, trigApprox, polyval_vector

   interface polyval
      module procedure polyval_scalar
      module procedure polyval_vector
      module procedure polyval_complex_scalar
   end interface polyval

contains

   real(dp) function polyval_scalar(p,x)
      real(dp),intent(in)::p(:),x
      integer::i
      if(size(p)==0)then
         polyval_scalar=0.0_dp; return
      end if
      polyval_scalar=p(1)
      do i=2,size(p)
         polyval_scalar=polyval_scalar*x+p(i)
      end do
   end function polyval_scalar

   function polyval_vector(p,x) result(y)
      real(dp),intent(in)::p(:),x(:)
      real(dp),allocatable::y(:)
      integer::i
      allocate(y(size(x))); y=0.0_dp
      if(size(p)==0)return
      y=p(1)
      do i=2,size(p)
         y=y*x+p(i)
      end do
   end function polyval_vector

   complex(dp) function polyval_complex_scalar(p,x)
      complex(dp),intent(in)::p(:),x
      integer::i
      if(size(p)==0)then
         polyval_complex_scalar=(0.0_dp,0.0_dp); return
      end if
      polyval_complex_scalar=p(1)
      do i=2,size(p)
         polyval_complex_scalar=polyval_complex_scalar*x+p(i)
      end do
   end function polyval_complex_scalar

   function polyvalm(p,a) result(y)
      real(dp),intent(in)::p(:),a(:,:)
      real(dp),allocatable::y(:,:),id(:,:)
      integer::n,i
      n=size(a,1); allocate(y(n,n),id(n,n)); id=0.0_dp
      do i=1,n; id(i,i)=1.0_dp; end do
      if(size(a,2)/=n .or. size(p)==0)then
         y=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      y=p(1)*id
      do i=2,size(p)
         y=matmul(y,a)+p(i)*id
      end do
   end function polyvalm

   function trim_leading(p,tolerance) result(q)
      real(dp),intent(in)::p(:)
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::q(:)
      real(dp)::tol
      integer::i
      tol=0.0_dp
      if(present(tolerance))tol=tolerance
      i=1
      do while(i<size(p) .and. abs(p(i))<=tol)
         i=i+1
      end do
      allocate(q(size(p)-i+1)); q=p(i:)
   end function trim_leading

   function polyadd(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp),allocatable::c(:)
      integer::n,oa,ob
      n=max(size(a),size(b)); allocate(c(n)); c=0.0_dp
      oa=n-size(a); ob=n-size(b)
      c(oa+1:)=c(oa+1:)+a
      c(ob+1:)=c(ob+1:)+b
      c=trim_leading(c)
   end function polyadd

   function polymul(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp),allocatable::c(:)
      integer::i,j
      if(size(a)==0 .or. size(b)==0)then
         allocate(c(0)); return
      end if
      allocate(c(size(a)+size(b)-1)); c=0.0_dp
      do i=1,size(a)
         do j=1,size(b)
            c(i+j-1)=c(i+j-1)+a(i)*b(j)
         end do
      end do
      c=trim_leading(c)
   end function polymul

   recursive function polypow(p,n) result(q)
      real(dp),intent(in)::p(:)
      integer,intent(in)::n
      real(dp),allocatable::q(:),h(:)
      if(n<0)then
         allocate(q(0))
      else if(n==0)then
         allocate(q(1)); q=1.0_dp
      else if(n==1)then
         allocate(q(size(p))); q=p
      else if(modulo(n,2)==0)then
         h=polypow(p,n/2); q=polymul(h,h)
      else
         h=polypow(p,n-1); q=polymul(h,p)
      end if
   end function polypow

   function polydiv(a,b) result(res)
      real(dp),intent(in)::a(:),b(:)
      type(polynomial_division_result)::res
      real(dp),allocatable::r(:),bb(:)
      integer::na,nb,i
      bb=trim_leading(b); r=trim_leading(a)
      na=size(r); nb=size(bb)
      if(nb==0 .or. abs(bb(1))<=tiny(1.0_dp))then
         allocate(res%quotient(0),res%remainder(size(r)))
         res%remainder=r; res%status=pracma_invalid_argument; return
      end if
      if(na<nb)then
         allocate(res%quotient(1),res%remainder(na))
         res%quotient=0.0_dp; res%remainder=r; return
      end if
      allocate(res%quotient(na-nb+1)); res%quotient=0.0_dp
      do i=1,na-nb+1
         res%quotient(i)=r(i)/bb(1)
         r(i:i+nb-1)=r(i:i+nb-1)-res%quotient(i)*bb
      end do
      if(nb>1)then
         allocate(res%remainder(nb-1)); res%remainder=r(na-nb+2:na)
         res%remainder=trim_leading(res%remainder,100.0_dp*eps_dp)
      else
         allocate(res%remainder(1)); res%remainder=0.0_dp
      end if
      res%status=pracma_ok
   end function polydiv

   function polyder(p) result(d)
      real(dp),intent(in)::p(:)
      real(dp),allocatable::d(:)
      integer::n,i
      n=size(p)-1
      if(n<=0)then
         allocate(d(1)); d=0.0_dp; return
      end if
      allocate(d(n))
      do i=1,n
         d(i)=real(n-i+1,dp)*p(i)
      end do
   end function polyder

   function polyint(p,c0) result(q)
      real(dp),intent(in)::p(:)
      real(dp),intent(in),optional::c0
      real(dp),allocatable::q(:)
      integer::n,i
      n=size(p); allocate(q(n+1))
      do i=1,n
         q(i)=p(i)/real(n-i+1,dp)
      end do
      q(n+1)=0.0_dp
      if(present(c0))q(n+1)=c0
   end function polyint

   function polytrans(p,shift,scale) result(q)
      real(dp),intent(in)::p(:)
      real(dp),intent(in),optional::shift,scale
      real(dp),allocatable::q(:),term(:),base(:)
      real(dp)::s,c
      integer::n,k
      s=0.0_dp; c=1.0_dp
      if(present(shift))s=shift
      if(present(scale))c=scale
      n=size(p)-1; allocate(q(1)); q=0.0_dp
      base=[c,s]
      do k=0,n
         term=polypow(base,n-k)
         term=p(k+1)*term
         q=polyadd(q,term)
      end do
   end function polytrans

   function poly2str(p,variable) result(text)
      real(dp),intent(in)::p(:)
      character(len=*),intent(in),optional::variable
      character(len=:),allocatable::text
      character(len=32)::buf
      character(len=16)::v
      integer::i,n,pow
      v='x'
      if(present(variable))v=variable
      text=''; n=size(p)-1
      do i=1,size(p)
         if(abs(p(i))<=100.0_dp*eps_dp)cycle
         write(buf,'(es16.8)')abs(p(i))
         if(len(text)>0)then
            if(p(i)>=0.0_dp)then; text=text//' + '; else; text=text//' - '; end if
         else if(p(i)<0.0_dp)then
            text='-'
         end if
         text=text//trim(adjustl(buf))
         pow=n-i+1
         if(pow>=1)text=text//'*'//trim(v)
         if(pow>=2)then
            write(buf,'(i0)')pow; text=text//'^'//trim(buf)
         end if
      end do
      if(len(text)==0)text='0'
   end function poly2str

   function roots(p,tolerance,max_iter,status) result(z)
      real(dp),intent(in)::p(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      integer,intent(out),optional::status
      complex(dp),allocatable::z(:)
      complex(dp),allocatable::coeff(:),znew(:)
      complex(dp)::prod,delta
      real(dp)::tol,radius,change
      integer::n,i,j,iter,niter,istat
      real(dp),allocatable::pp(:)
      pp=trim_leading(p,0.0_dp); n=size(pp)-1; istat=pracma_ok
      if(n<=0 .or. abs(pp(1))<=tiny(1.0_dp))then
         allocate(z(0)); istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      allocate(z(n),znew(n),coeff(n+1))
      coeff=cmplx(pp/pp(1),0.0_dp,dp)
      radius=1.0_dp+maxval(abs(pp(2:)/pp(1)))
      do i=1,n
         z(i)=radius*exp(cmplx(0.0_dp,2.0_dp*pi_dp*real(i-1,dp)/real(n,dp),dp))
      end do
      tol=1.0e-12_dp
      if(present(tolerance))tol=tolerance
      niter=2000
      if(present(max_iter))niter=max_iter
      do iter=1,niter
         change=0.0_dp
         do i=1,n
            prod=(1.0_dp,0.0_dp)
            do j=1,n
               if(j/=i)prod=prod*(z(i)-z(j))
            end do
            if(abs(prod)<=tiny(1.0_dp))prod=prod+cmplx(tol,tol,dp)
            delta=polyval_complex_scalar(coeff,z(i))/prod
            znew(i)=z(i)-delta
            change=max(change,abs(delta))
         end do
         z=znew
         if(change<=tol*(1.0_dp+maxval(abs(z))))exit
      end do
      if(iter>niter)istat=pracma_not_converged
      call sort_complex_roots(z)
      if(present(status))status=istat
   end function roots

   subroutine sort_complex_roots(z)
      complex(dp),intent(inout)::z(:)
      integer::i,j
      complex(dp)::v
      do i=2,size(z)
         v=z(i); j=i-1
         do while(j>=1)
            if(real(z(j),dp)<real(v,dp))exit
            if(real(z(j),dp)==real(v,dp) .and. aimag(z(j))<=aimag(v))exit
            z(j+1)=z(j); j=j-1
         end do
         z(j+1)=v
      end do
   end subroutine sort_complex_roots

   function rootsmult(p,tolerance) result(mult)
      real(dp),intent(in)::p(:)
      real(dp),intent(in),optional::tolerance
      integer,allocatable::mult(:)
      complex(dp),allocatable::z(:)
      real(dp)::tol
      integer::i,j
      tol=1.0e-7_dp
      if(present(tolerance))tol=tolerance
      z=roots(p); allocate(mult(size(z))); mult=1
      do i=1,size(z)
         do j=1,size(z)
            if(i/=j .and. abs(z(i)-z(j))<=tol*(1.0_dp+abs(z(i))))mult(i)=mult(i)+1
         end do
      end do
   end function rootsmult

   function polyroots(p) result(z)
      real(dp),intent(in)::p(:)
      complex(dp),allocatable::z(:)
      z=roots(p)
   end function polyroots

   function compan(p) result(a)
      real(dp),intent(in)::p(:)
      real(dp),allocatable::a(:,:)
      integer::n,i
      n=size(p)-1; allocate(a(max(0,n),max(0,n))); a=0.0_dp
      if(n<=0 .or. abs(p(1))<=tiny(1.0_dp))return
      a(1,:)=-p(2:)/p(1)
      do i=2,n
         a(i,i-1)=1.0_dp
      end do
   end function compan

   subroutine horner(p,x,value,derivative)
      real(dp),intent(in)::p(:),x
      real(dp),intent(out)::value
      real(dp),intent(out),optional::derivative
      real(dp)::d
      integer::i
      if(size(p)==0)then
         value=0.0_dp
         if(present(derivative))derivative=0.0_dp
         return
      end if
      value=p(1); d=0.0_dp
      do i=2,size(p)
         d=d*x+value; value=value*x+p(i)
      end do
      if(present(derivative))derivative=d
   end subroutine horner

   subroutine hornerdefl(p,x,value,q)
      real(dp),intent(in)::p(:),x
      real(dp),intent(out)::value
      real(dp),allocatable,intent(out)::q(:)
      integer::i,n
      n=size(p)-1
      if(n<=0)then
         allocate(q(0)); value=merge(p(1),0.0_dp,size(p)>0); return
      end if
      allocate(q(n)); q(1)=p(1)
      do i=2,n
         q(i)=q(i-1)*x+p(i)
      end do
      value=q(n)*x+p(n+1)
   end subroutine hornerdefl

   function polyfit(x,y,degree,status) result(res)
      real(dp),intent(in)::x(:),y(:)
      integer,intent(in)::degree
      integer,intent(out),optional::status
      type(polynomial_fit_result)::res
      real(dp),allocatable::v(:,:),coef(:)
      integer::i,j,istat
      if(size(x)/=size(y) .or. degree<0 .or. size(x)<degree+1)then
         allocate(res%coefficients(0)); res%status=pracma_invalid_argument
         if(present(status))status=res%status
         return
      end if
      allocate(v(size(x),degree+1))
      do i=1,size(x)
         do j=1,degree+1
            v(i,j)=x(i)**(degree-j+1)
         end do
      end do
      coef=mldivide(v,y,istat)
      allocate(res%coefficients(size(coef))); res%coefficients=coef
      res%residual_norm=sqrt(sum((matmul(v,coef)-y)**2))
      res%rank=degree+1; res%status=istat
      if(present(status))status=istat
   end function polyfit

   function polyfix(x,y,degree,fixed_powers,fixed_values,status) result(res)
      real(dp),intent(in)::x(:),y(:)
      integer,intent(in)::degree
      integer,intent(in)::fixed_powers(:)
      real(dp),intent(in)::fixed_values(:)
      integer,intent(out),optional::status
      type(polynomial_fit_result)::res
      logical,allocatable::is_fixed(:)
      integer,allocatable::free_powers(:)
      real(dp),allocatable::adjusted(:),v(:,:),coef_free(:)
      integer::i,j,k,istat
      if(size(fixed_powers)/=size(fixed_values) .or. size(x)/=size(y) .or. degree<0)then
         allocate(res%coefficients(0)); res%status=pracma_invalid_argument
         if(present(status))status=res%status
         return
      end if
      allocate(is_fixed(0:degree)); is_fixed=.false.
      do i=1,size(fixed_powers)
         if(fixed_powers(i)>=0 .and. fixed_powers(i)<=degree)is_fixed(fixed_powers(i))=.true.
      end do
      allocate(free_powers(count(.not.is_fixed))); k=0
      do i=0,degree
         if(.not.is_fixed(i))then; k=k+1; free_powers(k)=i; end if
      end do
      allocate(adjusted(size(y))); adjusted=y
      do i=1,size(fixed_powers)
         adjusted=adjusted-fixed_values(i)*x**fixed_powers(i)
      end do
      allocate(v(size(x),size(free_powers)))
      do j=1,size(free_powers)
         v(:,j)=x**free_powers(j)
      end do
      coef_free=mldivide(v,adjusted,istat)
      allocate(res%coefficients(degree+1)); res%coefficients=0.0_dp
      do j=1,size(free_powers)
         res%coefficients(degree-free_powers(j)+1)=coef_free(j)
      end do
      do i=1,size(fixed_powers)
         res%coefficients(degree-fixed_powers(i)+1)=fixed_values(i)
      end do
      res%residual_norm=sqrt(sum((polyval_vector(res%coefficients,x)-y)**2))
      res%rank=size(free_powers); res%status=istat
      if(present(status))status=istat
   end function polyfix

   function chebPoly(n,x) result(v)
      integer,intent(in)::n
      real(dp),intent(in)::x(:)
      real(dp),allocatable::v(:,:)
      integer::k
      allocate(v(size(x),n+1)); v(:,1)=1.0_dp
      if(n>=1)v(:,2)=x
      do k=2,n
         v(:,k+1)=2.0_dp*x*v(:,k)-v(:,k-1)
      end do
   end function chebPoly

   function chebCoeff(f,a,b,n) result(c)
      interface
         function f(x) result(y)
            import dp
            real(dp),intent(in)::x
            real(dp)::y
         end function f
      end interface
      real(dp),intent(in)::a,b
      integer,intent(in)::n
      real(dp),allocatable::c(:),theta(:),x(:),fx(:)
      integer::j,k
      allocate(c(n+1),theta(n+1),x(n+1),fx(n+1))
      do j=0,n
         theta(j+1)=pi_dp*(real(j,dp)+0.5_dp)/real(n+1,dp)
         x(j+1)=0.5_dp*(a+b)+0.5_dp*(b-a)*cos(theta(j+1))
         fx(j+1)=f(x(j+1))
      end do
      do k=0,n
         c(k+1)=2.0_dp/real(n+1,dp)*sum(fx*cos(real(k,dp)*theta))
      end do
      c(1)=0.5_dp*c(1)
   end function chebCoeff

   function chebApprox(c,a,b,x) result(y)
      real(dp),intent(in)::c(:),a,b,x(:)
      real(dp),allocatable::y(:),t(:),b1(:),b2(:),bk(:)
      integer::k
      allocate(y(size(x)),t(size(x)),b1(size(x)),b2(size(x)),bk(size(x)))
      t=(2.0_dp*x-a-b)/(b-a); b1=0.0_dp; b2=0.0_dp
      do k=size(c),2,-1
         bk=2.0_dp*t*b1-b2+c(k); b2=b1; b1=bk
      end do
      y=t*b1-b2+c(1)
   end function chebApprox

   function legendre(n,x) result(p)
      integer,intent(in)::n
      real(dp),intent(in)::x(:)
      real(dp),allocatable::p(:,:)
      integer::k
      allocate(p(size(x),n+1)); p(:,1)=1.0_dp
      if(n>=1)p(:,2)=x
      do k=2,n
         p(:,k+1)=((2.0_dp*real(k,dp)-1.0_dp)*x*p(:,k)-real(k-1,dp)*p(:,k-1))/real(k,dp)
      end do
   end function legendre

   function laguerre(n,x,alpha) result(p)
      integer,intent(in)::n
      real(dp),intent(in)::x(:)
      real(dp),intent(in),optional::alpha
      real(dp),allocatable::p(:,:)
      real(dp)::a
      integer::k
      a=0.0_dp
      if(present(alpha))a=alpha
      allocate(p(size(x),n+1)); p(:,1)=1.0_dp
      if(n>=1)p(:,2)=1.0_dp+a-x
      do k=2,n
         p(:,k+1)=((2.0_dp*real(k,dp)-1.0_dp+a-x)*p(:,k)- &
                    (real(k-1,dp)+a)*p(:,k-1))/real(k,dp)
      end do
   end function laguerre

   function bernstein(n,k,x) result(v)
      integer,intent(in)::n,k
      real(dp),intent(in)::x(:)
      real(dp),allocatable::v(:)
      allocate(v(size(x)))
      if(k<0 .or. k>n)then
         v=0.0_dp
      else
         v=binomial_real(n,k)*x**k*(1.0_dp-x)**(n-k)
      end if
   end function bernstein

   function bernsteinb(n,x) result(b)
      integer,intent(in)::n
      real(dp),intent(in)::x(:)
      real(dp),allocatable::b(:,:)
      integer::k
      allocate(b(size(x),n+1))
      do k=0,n
         b(:,k+1)=bernstein(n,k,x)
      end do
   end function bernsteinb

   pure real(dp) function binomial_real(n,k)
      integer,intent(in)::n,k
      integer::i,kk
      binomial_real=1.0_dp
      if(k<0 .or. k>n)then; binomial_real=0.0_dp; return; end if
      kk=min(k,n-k)
      do i=1,kk
         binomial_real=binomial_real*real(n-kk+i,dp)/real(i,dp)
      end do
   end function binomial_real

   subroutine pade(series,m,n,num,den,status)
      real(dp),intent(in)::series(:)
      integer,intent(in)::m,n
      real(dp),allocatable,intent(out)::num(:),den(:)
      integer,intent(out),optional::status
      real(dp),allocatable::a(:,:),rhs(:),q(:),p(:)
      integer::i,j,k,istat
      if(size(series)<m+n+1 .or. m<0 .or. n<0)then
         allocate(num(0),den(0)); istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      if(n==0)then
         allocate(den(1)); den=1.0_dp
         allocate(num(m+1)); num=series(1:m+1)
         if(present(status))status=pracma_ok
         return
      end if
      allocate(a(n,n),rhs(n))
      do i=1,n
         do j=1,n
            k=m+i-j
            if(k>=0)then; a(i,j)=series(k+1); else; a(i,j)=0.0_dp; end if
         end do
         rhs(i)=-series(m+i+1)
      end do
      q=mldivide(a,rhs,istat)
      allocate(den(n+1)); den(1)=1.0_dp; den(2:)=q
      allocate(p(m+1)); p=0.0_dp
      do i=0,m
         do j=0,min(i,n)
            p(i+1)=p(i+1)+den(j+1)*series(i-j+1)
         end do
      end do
      allocate(num(m+1)); num=p
      if(present(status))status=istat
   end subroutine pade

   function Poly(z) result(p)
      complex(dp),intent(in)::z(:)
      complex(dp),allocatable::p(:),q(:)
      integer::i
      allocate(p(1)); p=(1.0_dp,0.0_dp)
      do i=1,size(z)
         q=complex_polymul(p,[cmplx(1.0_dp,0.0_dp,dp),-z(i)])
         call move_alloc(q,p)
      end do
   end function Poly

   function complex_polymul(a,b) result(c)
      complex(dp),intent(in)::a(:),b(:)
      complex(dp),allocatable::c(:)
      integer::i,j
      allocate(c(size(a)+size(b)-1)); c=(0.0_dp,0.0_dp)
      do i=1,size(a); do j=1,size(b)
         c(i+j-1)=c(i+j-1)+a(i)*b(j)
      end do; end do
   end function complex_polymul

   function polygcf(p,q) result(g)
      real(dp),intent(in)::p(:),q(:)
      real(dp),allocatable::g(:),a(:),b(:)
      type(polynomial_division_result)::d
      real(dp)::scale
      a=trim_leading(p,1.0e-12_dp); b=trim_leading(q,1.0e-12_dp)
      do while(size(b)>1 .or. abs(b(1))>1.0e-10_dp)
         d=polydiv(a,b)
         if(size(d%remainder)==0)exit
         if(maxval(abs(d%remainder))<=1.0e-10_dp)exit
         a=b; b=d%remainder
      end do
      g=a
      if(size(g)>0)then
         scale=g(1)
         if(abs(scale)>tiny(1.0_dp))g=g/scale
      end if
   end function polygcf

   subroutine rationalfit(x,y,m,n,num,den,status)
      real(dp),intent(in)::x(:),y(:)
      integer,intent(in)::m,n
      real(dp),allocatable,intent(out)::num(:),den(:)
      integer,intent(out),optional::status
      real(dp),allocatable::a(:,:),rhs(:),coef(:)
      integer::i,j,istat
      if(size(x)/=size(y) .or. size(x)<m+n+1)then
         allocate(num(0),den(0)); istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      allocate(a(size(x),m+1+n),rhs(size(x))); rhs=y
      do i=1,size(x)
         do j=0,m
            a(i,j+1)=x(i)**(m-j)
         end do
         do j=1,n
            a(i,m+1+j)=-y(i)*x(i)**(n-j)
         end do
      end do
      coef=mldivide(a,rhs,istat)
      allocate(num(m+1),den(n+1))
      num=coef(1:m+1); den(1)=1.0_dp; den(2:)=coef(m+2:)
      if(present(status))status=istat
   end subroutine rationalfit

   function trigPoly(coefficients,x) result(y)
      real(dp),intent(in)::coefficients(:),x(:)
      real(dp),allocatable::y(:)
      integer::m,k
      m=(size(coefficients)-1)/2
      allocate(y(size(x))); y=coefficients(1)
      do k=1,m
         y=y+coefficients(2*k)*cos(real(k,dp)*x)+coefficients(2*k+1)*sin(real(k,dp)*x)
      end do
   end function trigPoly

   function trigApprox(x,y,degree,status) result(coefficients)
      real(dp),intent(in)::x(:),y(:)
      integer,intent(in)::degree
      integer,intent(out),optional::status
      real(dp),allocatable::coefficients(:),a(:,:)
      integer::k,istat
      if(size(x)/=size(y) .or. degree<0)then
         allocate(coefficients(0)); istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      allocate(a(size(x),2*degree+1)); a(:,1)=1.0_dp
      do k=1,degree
         a(:,2*k)=cos(real(k,dp)*x)
         a(:,2*k+1)=sin(real(k,dp)*x)
      end do
      coefficients=mldivide(a,y,istat)
      if(present(status))status=istat
   end function trigApprox

end module pracma_polynomial
