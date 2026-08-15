module ld_numerics
use ld_kinds, only: dp
use ld_interfaces, only: log_target_iface, vector_func_iface
implicit none
private
public :: numerical_gradient, numerical_hessian, numerical_jacobian, logsumexp, logadd
contains
function logsumexp(x) result(v)
   real(dp),intent(in)::x(:); real(dp)::v,m
   m=maxval(x)
   if(.not.(m>-huge(1.0_dp))) then; v=m; else; v=m+log(sum(exp(x-m))); end if
end function logsumexp
function logadd(x) result(v)
   real(dp),intent(in)::x(:); real(dp)::v
   v=logsumexp(x)
end function logadd
subroutine numerical_gradient(f,x,g,h)
   procedure(log_target_iface) :: f
   real(dp),intent(in)::x(:)
   real(dp),intent(out)::g(:)
   real(dp),intent(in),optional::h
   real(dp)::xp(size(x)),xm(size(x)),step
   integer::j
   step=1e-6_dp; if(present(h)) step=h
   do j=1,size(x)
      xp=x; xm=x; xp(j)=xp(j)+step; xm(j)=xm(j)-step
      g(j)=(f(xp)-f(xm))/(2.0_dp*step)
   end do
end subroutine numerical_gradient
subroutine numerical_hessian(f,x,hess,h)
   procedure(log_target_iface) :: f
   real(dp),intent(in)::x(:)
   real(dp),intent(out)::hess(:,:)
   real(dp),intent(in),optional::h
   real(dp)::xp(size(x)),xm(size(x)),xpp(size(x)),xpm(size(x)),xmp(size(x)),xmm(size(x)),step,f0
   integer::i,j
   step=1e-4_dp; if(present(h)) step=h; f0=f(x); hess=0.0_dp
   do i=1,size(x)
      xp=x; xm=x; xp(i)=xp(i)+step; xm(i)=xm(i)-step
      hess(i,i)=(f(xp)-2.0_dp*f0+f(xm))/(step*step)
      do j=i+1,size(x)
         xpp=x; xpm=x; xmp=x; xmm=x
         xpp(i)=xpp(i)+step; xpp(j)=xpp(j)+step
         xpm(i)=xpm(i)+step; xpm(j)=xpm(j)-step
         xmp(i)=xmp(i)-step; xmp(j)=xmp(j)+step
         xmm(i)=xmm(i)-step; xmm(j)=xmm(j)-step
         hess(i,j)=(f(xpp)-f(xpm)-f(xmp)+f(xmm))/(4.0_dp*step*step)
         hess(j,i)=hess(i,j)
      end do
   end do
end subroutine numerical_hessian
subroutine numerical_jacobian(f,x,m,jac,h)
   procedure(vector_func_iface) :: f
   real(dp),intent(in)::x(:)
   integer,intent(in)::m
   real(dp),intent(out)::jac(:,:)
   real(dp),intent(in),optional::h
   real(dp)::xp(size(x)),xm(size(x)),yp(m),ym(m),step
   integer::j
   step=1e-6_dp; if(present(h)) step=h
   do j=1,size(x)
      xp=x; xm=x; xp(j)=xp(j)+step; xm(j)=xm(j)-step
      call f(xp,yp); call f(xm,ym); jac(:,j)=(yp-ym)/(2.0_dp*step)
   end do
end subroutine numerical_jacobian
end module ld_numerics
