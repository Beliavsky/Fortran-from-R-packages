! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_differentiation
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use pracma_kinds, only : dp, eps_dp
   use pracma_status
   use pracma_callbacks
   implicit none
   private

   public :: gradient, grad, jacobian, hessian, hessvec, hessdiag, laplacian
   public :: numderiv, numdiff, complexstep, grad_csd, jacobian_csd
   public :: hessian_csd, laplacian_csd, fornberg, fderiv, taylor

contains

   function gradient(y, x) result(g)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in), optional :: x(:)
      real(dp), allocatable :: g(:)
      real(dp), allocatable :: xx(:)
      integer :: n, i
      n = size(y); allocate(g(n), xx(n))
      if (n == 0) return
      if (present(x)) then
         if (size(x) /= n) then
            g = ieee_value(0.0_dp, ieee_quiet_nan); return
         end if
         xx = x
      else
         xx = [(real(i,dp), i=1,n)]
      end if
      if (n == 1) then
         g = 0.0_dp
      else
         g(1) = (y(2)-y(1))/(xx(2)-xx(1))
         g(n) = (y(n)-y(n-1))/(xx(n)-xx(n-1))
         do i=2,n-1
            g(i) = ((xx(i)-xx(i-1))**2*(y(i+1)-y(i)) + &
                    (xx(i+1)-xx(i))**2*(y(i)-y(i-1))) / &
                   ((xx(i+1)-xx(i))*(xx(i)-xx(i-1))*(xx(i+1)-xx(i-1)))
         end do
      end if
   end function gradient

   subroutine grad(f, x, g, step, method, status)
      procedure(objective_function) :: f
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      real(dp), intent(in), optional :: step
      character(len=*), intent(in), optional :: method
      integer, intent(out), optional :: status
      real(dp), allocatable :: xp(:), xm(:)
      real(dp) :: h, f0, fp, fm
      integer :: i, istat
      character(len=16) :: meth
      istat=pracma_ok; meth='central'
      if(present(method))meth=trim(method)
      if(size(g)/=size(x))then
         g=0.0_dp; istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      allocate(xp(size(x)),xm(size(x))); xp=x; xm=x
      f0=f(x)
      do i=1,size(x)
         h=sqrt(eps_dp)*max(1.0_dp,abs(x(i)))
         if(present(step))h=step
         select case(meth)
         case('forward','simple')
            xp=x; xp(i)=x(i)+h; fp=f(xp); g(i)=(fp-f0)/h
         case('backward')
            xm=x; xm(i)=x(i)-h; fm=f(xm); g(i)=(f0-fm)/h
         case default
            h=eps_dp**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x(i)))
            if(present(step))h=step
            xp=x; xm=x; xp(i)=x(i)+h; xm(i)=x(i)-h
            fp=f(xp); fm=f(xm); g(i)=(fp-fm)/(2.0_dp*h)
         end select
         if(.not.ieee_is_finite(g(i)))istat=pracma_nonfinite
      end do
      if(present(status))status=istat
   end subroutine grad

   subroutine jacobian(f,x,m,jac,step,method,status)
      procedure(vector_function) :: f
      real(dp),intent(in)::x(:)
      integer,intent(in)::m
      real(dp),intent(out)::jac(:,:)
      real(dp),intent(in),optional::step
      character(len=*),intent(in),optional::method
      integer,intent(out),optional::status
      real(dp),allocatable::xp(:),xm(:),f0(:),fp(:),fm(:)
      real(dp)::h
      integer::i,istat
      character(len=16)::meth
      meth='central'; if(present(method))meth=trim(method)
      istat=pracma_ok
      if(any(shape(jac)/=[m,size(x)]))then
         jac=0.0_dp; istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      allocate(xp(size(x)),xm(size(x)),f0(m),fp(m),fm(m))
      call f(x,f0)
      do i=1,size(x)
         if(meth=='forward' .or. meth=='simple')then
            h=sqrt(eps_dp)*max(1.0_dp,abs(x(i)))
            if(present(step))h=step
            xp=x; xp(i)=x(i)+h; call f(xp,fp); jac(:,i)=(fp-f0)/h
         else
            h=eps_dp**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x(i)))
            if(present(step))h=step
            xp=x; xm=x; xp(i)=x(i)+h; xm(i)=x(i)-h
            call f(xp,fp); call f(xm,fm); jac(:,i)=(fp-fm)/(2.0_dp*h)
         end if
      end do
      if(any(.not.ieee_is_finite(jac)))istat=pracma_nonfinite
      if(present(status))status=istat
   end subroutine jacobian

   subroutine hessian(f,x,hess,step,status)
      procedure(objective_function)::f
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::hess(:,:)
      real(dp),intent(in),optional::step
      integer,intent(out),optional::status
      real(dp),allocatable::xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:)
      real(dp)::hi,hj,f0
      integer::i,j,n,istat
      n=size(x); istat=pracma_ok
      if(any(shape(hess)/=[n,n]))then
         hess=0.0_dp; istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      allocate(xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n)); f0=f(x)
      do i=1,n
         hi=eps_dp**0.25_dp*max(1.0_dp,abs(x(i)))
         if(present(step))hi=step
         xp=x; xm=x; xp(i)=x(i)+hi; xm(i)=x(i)-hi
         hess(i,i)=(f(xp)-2.0_dp*f0+f(xm))/(hi*hi)
         do j=i+1,n
            hj=eps_dp**0.25_dp*max(1.0_dp,abs(x(j)))
            if(present(step))hj=step
            xpp=x; xpm=x; xmp=x; xmm=x
            xpp(i)=x(i)+hi; xpp(j)=x(j)+hj
            xpm(i)=x(i)+hi; xpm(j)=x(j)-hj
            xmp(i)=x(i)-hi; xmp(j)=x(j)+hj
            xmm(i)=x(i)-hi; xmm(j)=x(j)-hj
            hess(i,j)=(f(xpp)-f(xpm)-f(xmp)+f(xmm))/(4.0_dp*hi*hj)
            hess(j,i)=hess(i,j)
         end do
      end do
      if(any(.not.ieee_is_finite(hess)))istat=pracma_nonfinite
      if(present(status))status=istat
   end subroutine hessian

   subroutine hessvec(f,x,v,hv,step,status)
      procedure(objective_function)::f
      real(dp),intent(in)::x(:),v(:)
      real(dp),intent(out)::hv(:)
      real(dp),intent(in),optional::step
      integer,intent(out),optional::status
      real(dp),allocatable::gp(:),gm(:),xp(:),xm(:)
      real(dp)::h
      integer::istat1,istat2,istat
      if(size(v)/=size(x) .or. size(hv)/=size(x))then
         hv=0.0_dp; istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      h=eps_dp**(1.0_dp/3.0_dp)/max(1.0_dp,sqrt(sum(v*v)))
      if(present(step))h=step
      allocate(gp(size(x)),gm(size(x)),xp(size(x)),xm(size(x)))
      xp=x+h*v; xm=x-h*v
      call grad(f,xp,gp,status=istat1); call grad(f,xm,gm,status=istat2)
      hv=(gp-gm)/(2.0_dp*h)
      istat=max(istat1,istat2)
      if(present(status))status=istat
   end subroutine hessvec

   subroutine hessdiag(f,x,d,step,status)
      procedure(objective_function)::f
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::d(:)
      real(dp),intent(in),optional::step
      integer,intent(out),optional::status
      real(dp),allocatable::h(:,:)
      integer::i,istat
      allocate(h(size(x),size(x))); call hessian(f,x,h,step,istat)
      do i=1,size(x); d(i)=h(i,i); end do
      if(present(status))status=istat
   end subroutine hessdiag

   function laplacian(f,x,step,status) result(v)
      procedure(objective_function)::f
      real(dp),intent(in)::x(:)
      real(dp),intent(in),optional::step
      integer,intent(out),optional::status
      real(dp)::v
      real(dp),allocatable::d(:)
      integer::istat
      allocate(d(size(x))); call hessdiag(f,x,d,step,istat); v=sum(d)
      if(present(status))status=istat
   end function laplacian

   function numderiv(f,x,order,step,status) result(v)
      procedure(scalar_function)::f
      real(dp),intent(in)::x
      integer,intent(in),optional::order
      real(dp),intent(in),optional::step
      integer,intent(out),optional::status
      real(dp)::v,h
      integer::n,istat
      n=1; if(present(order))n=order
      h=eps_dp**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x)); if(present(step))h=step
      select case(n)
      case(0); v=f(x)
      case(1); v=(f(x+h)-f(x-h))/(2.0_dp*h)
      case(2); v=(f(x+h)-2.0_dp*f(x)+f(x-h))/(h*h)
      case(3); v=(f(x+2*h)-2*f(x+h)+2*f(x-h)-f(x-2*h))/(2*h**3)
      case(4); v=(f(x-2*h)-4*f(x-h)+6*f(x)-4*f(x+h)+f(x+2*h))/h**4
      case default
         v=fderiv(f,x,n,h,istat)
         if(present(status))status=istat
         return
      end select
      istat=merge(pracma_ok,pracma_nonfinite,ieee_is_finite(v))
      if(present(status))status=istat
   end function numderiv

   function numdiff(f,x,order,step,status) result(v)
      procedure(scalar_function)::f
      real(dp),intent(in)::x
      integer,intent(in),optional::order
      real(dp),intent(in),optional::step
      integer,intent(out),optional::status
      real(dp)::v
      v=numderiv(f,x,order,step,status)
   end function numdiff

   function complexstep(f,x,h) result(v)
      procedure(complex_scalar_function)::f
      real(dp),intent(in)::x
      real(dp),intent(in),optional::h
      real(dp)::v,hh
      hh=1.0e-20_dp; if(present(h))hh=h
      v=aimag(f(cmplx(x,hh,dp)))/hh
   end function complexstep

   subroutine grad_csd(f,x,g,h,status)
      procedure(complex_objective_function)::f
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::g(:)
      real(dp),intent(in),optional::h
      integer,intent(out),optional::status
      complex(dp),allocatable::z(:)
      real(dp)::hh
      integer::i,istat
      hh=1.0e-20_dp; if(present(h))hh=h
      if(size(g)/=size(x))then
         g=0.0_dp; istat=pracma_dimension_mismatch
      else
         allocate(z(size(x))); z=cmplx(x,0.0_dp,dp)
         do i=1,size(x)
            z=cmplx(x,0.0_dp,dp); z(i)=cmplx(x(i),hh,dp)
            g(i)=aimag(f(z))/hh
         end do
         istat=pracma_ok
      end if
      if(present(status))status=istat
   end subroutine grad_csd

   subroutine jacobian_csd(f,x,m,jac,h,status)
      procedure(complex_vector_function)::f
      real(dp),intent(in)::x(:)
      integer,intent(in)::m
      real(dp),intent(out)::jac(:,:)
      real(dp),intent(in),optional::h
      integer,intent(out),optional::status
      complex(dp),allocatable::z(:),y(:)
      real(dp)::hh
      integer::i,istat
      hh=1.0e-20_dp; if(present(h))hh=h
      if(any(shape(jac)/=[m,size(x)]))then
         jac=0.0_dp; istat=pracma_dimension_mismatch
      else
         allocate(z(size(x)),y(m))
         do i=1,size(x)
            z=cmplx(x,0.0_dp,dp); z(i)=cmplx(x(i),hh,dp)
            call f(z,y); jac(:,i)=aimag(y)/hh
         end do
         istat=pracma_ok
      end if
      if(present(status))status=istat
   end subroutine jacobian_csd

   subroutine hessian_csd(f,x,hess,h,status)
      procedure(complex_objective_function)::f
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::hess(:,:)
      real(dp),intent(in),optional::h
      integer,intent(out),optional::status
      real(dp)::hh,f0
      complex(dp),allocatable::z(:)
      integer::i,j,n,istat
      n=size(x); hh=1.0e-5_dp; if(present(h))hh=h
      if(any(shape(hess)/=[n,n]))then
         hess=0.0_dp; istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      allocate(z(n)); f0=real(f(cmplx(x,0.0_dp,dp)),dp)
      do i=1,n
         z=cmplx(x,0.0_dp,dp); z(i)=cmplx(x(i),hh,dp)
         hess(i,i)=2.0_dp*(f0-real(f(z),dp))/(hh*hh)
         do j=i+1,n
            z=cmplx(x,0.0_dp,dp); z(i)=cmplx(x(i),hh,dp); z(j)=cmplx(x(j),hh,dp)
            hess(i,j)=-(real(f(z),dp)-f0+0.5_dp*hh*hh*(hess(i,i)+hess(j,j)))/(hh*hh)
            hess(j,i)=hess(i,j)
         end do
      end do
      istat=pracma_ok
      if(present(status))status=istat
   end subroutine hessian_csd

   function laplacian_csd(f,x,h,status) result(v)
      procedure(complex_objective_function)::f
      real(dp),intent(in)::x(:)
      real(dp),intent(in),optional::h
      integer,intent(out),optional::status
      real(dp)::v
      real(dp),allocatable::hs(:,:)
      integer::i,istat
      allocate(hs(size(x),size(x))); call hessian_csd(f,x,hs,h,istat)
      v=0.0_dp; do i=1,size(x); v=v+hs(i,i); end do
      if(present(status))status=istat
   end function laplacian_csd

   subroutine fornberg(x0,x,m,c,status)
      real(dp),intent(in)::x0,x(:)
      integer,intent(in)::m
      real(dp),intent(out)::c(:,:)
      integer,intent(out),optional::status
      real(dp)::c1,c2,c3,c4,c5
      integer::n,i,j,k,mn,istat
      n=size(x)
      if(size(c,1)/=n .or. size(c,2)/=m+1 .or. m<0 .or. m>=n)then
         c=0.0_dp; istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      c=0.0_dp; c(1,1)=1.0_dp; c1=1.0_dp; c4=x(1)-x0
      do i=2,n
         mn=min(i-1,m); c2=1.0_dp; c5=c4; c4=x(i)-x0
         do j=1,i-1
            c3=x(i)-x(j); c2=c2*c3
            if(j==i-1)then
               do k=mn,1,-1
                  c(i,k+1)=c1*(real(k,dp)*c(i-1,k)-c5*c(i-1,k+1))/c2
               end do
               c(i,1)=-c1*c5*c(i-1,1)/c2
            end if
            do k=mn,1,-1
               c(j,k+1)=(c4*c(j,k+1)-real(k,dp)*c(j,k))/c3
            end do
            c(j,1)=c4*c(j,1)/c3
         end do
         c1=c2
      end do
      istat=pracma_ok
      if(present(status))status=istat
   end subroutine fornberg

   function fderiv(f,x,order,step,status) result(v)
      procedure(scalar_function)::f
      real(dp),intent(in)::x
      integer,intent(in)::order
      real(dp),intent(in),optional::step
      integer,intent(out),optional::status
      real(dp)::v,h
      real(dp),allocatable::nodes(:),c(:,:),values(:)
      integer::n,i,istat
      if(order<0)then
         v=ieee_value(0.0_dp,ieee_quiet_nan); istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      n=max(order+3,7)
      if(modulo(n,2)==0)n=n+1
      h=eps_dp**(1.0_dp/real(order+2,dp))*max(1.0_dp,abs(x))
      if(present(step))h=step
      allocate(nodes(n),c(n,order+1),values(n))
      do i=1,n
         nodes(i)=x+real(i-(n+1)/2,dp)*h; values(i)=f(nodes(i))
      end do
      call fornberg(x,nodes,order,c,istat); v=dot_product(c(:,order+1),values)
      if(present(status))status=istat
   end function fderiv

   function taylor(f,x0,degree,step,status) result(coefficients)
      procedure(scalar_function)::f
      real(dp),intent(in)::x0
      integer,intent(in)::degree
      real(dp),intent(in),optional::step
      integer,intent(out),optional::status
      real(dp),allocatable::coefficients(:)
      real(dp)::fact
      integer::k,istat,current
      if(degree<0)then
         allocate(coefficients(0)); current=pracma_invalid_argument
      else
         allocate(coefficients(degree+1)); fact=1.0_dp; current=pracma_ok
         coefficients(1)=f(x0)
         do k=1,degree
            fact=fact*real(k,dp)
            coefficients(k+1)=fderiv(f,x0,k,step,istat)/fact
            if(istat/=pracma_ok)current=istat
         end do
      end if
      if(present(status))status=current
   end function taylor

end module pracma_differentiation
