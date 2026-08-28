module lavaan_optimizer
   use lavaan_kinds, only : dp
   implicit none
   private
   public :: bfgs_minimize
   abstract interface
      function scalar_fun(x) result(f)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: f
      end function scalar_fun
   end interface
contains
   subroutine numeric_gradient(fun,x,g)
      procedure(scalar_fun)::fun
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::g(:)
      real(dp),allocatable::xp(:),xm(:)
      real(dp)::h
      integer::i
      allocate(xp(size(x)),xm(size(x)))
      do i=1,size(x)
         h=1.0e-5_dp*max(1.0_dp,abs(x(i)))
         xp=x
         xm=x
         xp(i)=xp(i)+h
         xm(i)=xm(i)-h
         g(i)=(fun(xp)-fun(xm))/(2*h)
      end do
   end subroutine numeric_gradient

   subroutine bfgs_minimize(fun,x,fval,converged,iterations,maxiter,tol)
      procedure(scalar_fun)::fun
      real(dp),intent(inout)::x(:)
      real(dp),intent(out)::fval
      logical,intent(out)::converged
      integer,intent(out)::iterations
      integer,intent(in),optional::maxiter
      real(dp),intent(in),optional::tol
      integer::n,i,it,mx
      real(dp)::ft,fn,alpha,gtol,ys,rho
      real(dp),allocatable::h(:,:),g(:),gn(:),p(:),xn(:),s(:),y(:),eye(:,:),hy(:)
      mx=500
      if(present(maxiter)) mx=maxiter
      gtol=1e-7_dp
      if(present(tol)) gtol=tol
      n=size(x)
      allocate(h(n,n),g(n),gn(n),p(n),xn(n),s(n),y(n),eye(n,n),hy(n))
      h=0
      eye=0
      do i=1,n
      h(i,i)=1
      eye(i,i)=1
      end do
      ft=fun(x)
      call numeric_gradient(fun,x,g)
      converged=.false.
      do it=1,mx
         if(maxval(abs(g))<gtol) then
         converged=.true.
         exit
         end if
         p=-matmul(h,g)
         if(dot_product(p,g)>=0) p=-g
         alpha=1.0_dp
         do
            xn=x+alpha*p
            fn=fun(xn)
            if(fn < ft + 1.0e-4_dp*alpha*dot_product(g,p)) exit
            alpha=alpha*0.5_dp
            if(alpha<1.0e-10_dp) exit
         end do
         if(alpha<1.0e-10_dp) exit
         call numeric_gradient(fun,xn,gn)
         s=xn-x
         y=gn-g
         ys=dot_product(y,s)
         if(ys>1.0e-12_dp) then
            rho=1.0_dp/ys
            h=matmul(eye-rho*spread(s,2,n)*spread(y,1,n), &
               matmul(h,eye-rho*spread(y,2,n)*spread(s,1,n))) + rho*spread(s,2,n)*spread(s,1,n)
         else
            h=eye
         end if
         x=xn
         g=gn
         ft=fn
      end do
      iterations=min(it,mx)
      fval=ft
   end subroutine bfgs_minimize
end module lavaan_optimizer
