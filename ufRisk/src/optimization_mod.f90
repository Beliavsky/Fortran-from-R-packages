! SPDX-License-Identifier: GPL-3.0-only
module optimization_mod
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use kind_mod, only: dp
   implicit none
   private
   type, public :: optimization_result_t
      real(dp), allocatable :: parameters(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type optimization_result_t
   abstract interface
      pure function objective_interface(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_interface
   end interface
   public :: bfgs_minimize_fd, finite_difference_hessian
contains
   pure function bfgs_minimize_fd(objective, start, max_iterations, tolerance) result(out)
      procedure(objective_interface) :: objective
      real(dp), intent(in) :: start(:)
      integer, intent(in) :: max_iterations
      real(dp), intent(in) :: tolerance
      type(optimization_result_t) :: out
      real(dp), allocatable :: x(:),g(:),gnew(:),h(:,:),direction(:),trial(:),s(:),y(:),hy(:)
      real(dp) :: fx,ftrial,step,ys,yhy,rho
      integer :: n,i,iter,ls
      n=size(start); allocate(x(n),g(n),gnew(n),h(n,n),direction(n),trial(n),s(n),y(n),hy(n))
      x=start; fx=objective(x); call gradient_fd(objective,x,g)
      h=0.0_dp;do i=1,n;h(i,i)=1.0_dp;end do
      do iter=1,max_iterations
         if(maxval(abs(g))<=tolerance*(1.0_dp+abs(fx)))then
            out%converged=.true.;exit
         end if
         direction=-matmul(h,g)
         if(dot_product(direction,g)>=0.0_dp .or. any(.not.(abs(direction)<=huge(1.0_dp))))direction=-g
         step=1.0_dp
         do ls=1,40
            trial=x+step*direction;ftrial=objective(trial)
            if(ftrial<huge(1.0_dp)/10.0_dp .and. ftrial<=fx+1.0e-4_dp*step*dot_product(g,direction))exit
            step=0.5_dp*step
         end do
         if(ls>40)then;out%info=2;exit;end if
         call gradient_fd(objective,trial,gnew)
         s=trial-x;y=gnew-g;ys=dot_product(y,s)
         if(ys>sqrt(epsilon(1.0_dp))*sqrt(max(dot_product(y,y)*dot_product(s,s),tiny(1.0_dp))))then
            hy=matmul(h,y);yhy=dot_product(y,hy);rho=1.0_dp/ys
            h=h + ((ys+yhy)*rho*rho)*outer(s,s)-rho*(outer(hy,s)+outer(s,hy))
         else
            h=0.0_dp;do i=1,n;h(i,i)=1.0_dp;end do
         end if
         x=trial;fx=ftrial;g=gnew
      end do
      out%parameters=x;out%value=fx;out%iterations=min(iter,max_iterations)
      if(.not.out%converged .and. out%info==0)out%info=1
   contains
      pure function outer(a,b) result(c)
         real(dp),intent(in)::a(:),b(:);real(dp)::c(size(a),size(b));integer::ii,jj
         do jj=1,size(b);do ii=1,size(a);c(ii,jj)=a(ii)*b(jj);end do;end do
      end function outer
   end function bfgs_minimize_fd

   pure subroutine gradient_fd(objective,x,g)
      procedure(objective_interface)::objective
      real(dp),intent(in)::x(:);real(dp),intent(out)::g(:)
      real(dp)::xp(size(x)),xm(size(x)),step,fp,fm,f0
      integer::i
      do i=1,size(x)
         step=epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x(i)))
         xp=x;xm=x;xp(i)=xp(i)+step;xm(i)=xm(i)-step
         fp=objective(xp);fm=objective(xm)
         if(ieee_is_finite(fp).and.ieee_is_finite(fm).and. &
            abs(fp)<huge(1.0_dp)/10.0_dp.and.abs(fm)<huge(1.0_dp)/10.0_dp)then
            g(i)=(fp-fm)/(2.0_dp*step)
         else
            f0=objective(x)
            if(ieee_is_finite(fp).and.abs(fp)<huge(1.0_dp)/10.0_dp.and. &
               ieee_is_finite(f0).and.abs(f0)<huge(1.0_dp)/10.0_dp)then
               g(i)=(fp-f0)/step
            else if(ieee_is_finite(fm).and.abs(fm)<huge(1.0_dp)/10.0_dp.and. &
               ieee_is_finite(f0).and.abs(f0)<huge(1.0_dp)/10.0_dp)then
               g(i)=(f0-fm)/step
            else
               g(i)=0.0_dp
            end if
         end if
      end do
   end subroutine gradient_fd

   pure function finite_difference_hessian(objective,x) result(hessian)
      procedure(objective_interface)::objective
      real(dp),intent(in)::x(:)
      real(dp),allocatable::hessian(:,:)
      real(dp)::xpp(size(x)),xpm(size(x)),xmp(size(x)),xmm(size(x)),hi,hj,f0,fpp,fpm,fmp,fmm
      integer::i,j,n
      n=size(x);allocate(hessian(n,n));hessian=0.0_dp;f0=objective(x)
      do i=1,n
         hi=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(i)))
         xpp=x;xmm=x;xpp(i)=xpp(i)+hi;xmm(i)=xmm(i)-hi
         fpp=objective(xpp);fmm=objective(xmm)
         if(ieee_is_finite(fpp).and.ieee_is_finite(fmm).and.ieee_is_finite(f0).and. &
            max(abs(fpp),abs(fmm),abs(f0))<huge(1.0_dp)/10.0_dp)then
            hessian(i,i)=(fpp-2.0_dp*f0+fmm)/(hi*hi)
         else
            hessian(i,i)=0.0_dp
         end if
         do j=i+1,n
            hj=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(j)))
            xpp=x;xpm=x;xmp=x;xmm=x
            xpp(i)=xpp(i)+hi;xpp(j)=xpp(j)+hj
            xpm(i)=xpm(i)+hi;xpm(j)=xpm(j)-hj
            xmp(i)=xmp(i)-hi;xmp(j)=xmp(j)+hj
            xmm(i)=xmm(i)-hi;xmm(j)=xmm(j)-hj
            fpp=objective(xpp);fpm=objective(xpm);fmp=objective(xmp);fmm=objective(xmm)
            if(ieee_is_finite(fpp).and.ieee_is_finite(fpm).and.ieee_is_finite(fmp).and. &
               ieee_is_finite(fmm).and.max(abs(fpp),abs(fpm),abs(fmp),abs(fmm))<huge(1.0_dp)/10.0_dp)then
               hessian(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj)
            else
               hessian(i,j)=0.0_dp
            end if
            hessian(j,i)=hessian(i,j)
         end do
      end do
   end function finite_difference_hessian
end module optimization_mod
