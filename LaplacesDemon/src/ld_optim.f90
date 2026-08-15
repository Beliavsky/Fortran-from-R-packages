module ld_optim
use ld_kinds, only: dp
use ld_interfaces, only: log_target_iface
use ld_linalg, only: inverse_spd, make_positive_definite, outer_product
use ld_numerics, only: numerical_gradient, numerical_hessian
implicit none
private
public :: optim_result_t, laplace_result_t, bfgs_maximize, laplace_approximation

type :: optim_result_t
   real(dp), allocatable :: par(:)
   real(dp) :: value=-huge(1.0_dp)
   integer :: iterations=0
   logical :: converged=.false.
end type

type :: laplace_result_t
   real(dp), allocatable :: mode(:), covariance(:,:), hessian(:,:)
   real(dp) :: log_posterior=-huge(1.0_dp)
   integer :: iterations=0
   logical :: converged=.false.
end type
contains
subroutine bfgs_maximize(f,x0,res,max_iter,tol)
   procedure(log_target_iface) :: f
   real(dp),intent(in)::x0(:)
   type(optim_result_t),intent(out)::res
   integer,intent(in),optional::max_iter
   real(dp),intent(in),optional::tol
   integer::n,it,maxit,ls
   real(dp)::eps,alpha,fx,fn,ys,rho
   real(dp),allocatable::x(:),xn(:),g(:),gn(:),p(:),s(:),y(:),h(:,:),iunit(:,:)
   n=size(x0); maxit=500; if(present(max_iter)) maxit=max_iter; eps=1e-7_dp; if(present(tol)) eps=tol
   allocate(x(n),xn(n),g(n),gn(n),p(n),s(n),y(n),h(n,n),iunit(n,n),res%par(n))
   x=x0; h=0.0_dp; iunit=0.0_dp
   do it=1,n; h(it,it)=1.0_dp; iunit(it,it)=1.0_dp; end do
   fx=f(x); call numerical_gradient(f,x,g)
   do it=1,maxit
      if(maxval(abs(g))<eps) exit
      p=matmul(h,g)
      if(dot_product(g,p)<=0.0_dp) p=g
      alpha=1.0_dp
      do ls=1,30
         xn=x+alpha*p; fn=f(xn)
         if(fn>=fx+1e-4_dp*alpha*dot_product(g,p)) exit
         alpha=0.5_dp*alpha
      end do
      if(alpha<1e-12_dp) exit
      call numerical_gradient(f,xn,gn)
      s=xn-x; y=gn-g; ys=dot_product(y,s)
      if(abs(ys)>1e-14_dp) then
         ! Maximization BFGS with y = grad_new-grad_old: update inverse of -Hessian using -y.
         y=-y; ys=dot_product(y,s)
         if(ys>1e-14_dp) then
            rho=1.0_dp/ys
            h=matmul(iunit-rho*outer_product(s,y),matmul(h,iunit-rho*outer_product(y,s)))+rho*outer_product(s,s)
         end if
      end if
      x=xn; g=gn; fx=fn
      if(maxval(abs(s))<eps*(1.0_dp+maxval(abs(x)))) exit
   end do
   res%par=x; res%value=fx; res%iterations=min(it,maxit); res%converged=(maxval(abs(g))<sqrt(eps))
end subroutine bfgs_maximize

subroutine laplace_approximation(f,x0,res,max_iter,tol)
   procedure(log_target_iface) :: f
   real(dp),intent(in)::x0(:)
   type(laplace_result_t),intent(out)::res
   integer,intent(in),optional::max_iter
   real(dp),intent(in),optional::tol
   type(optim_result_t)::opt
   real(dp),allocatable::neg_h(:,:)
   integer::n,info
   call bfgs_maximize(f,x0,opt,max_iter,tol); n=size(x0)
   allocate(res%mode(n),res%covariance(n,n),res%hessian(n,n),neg_h(n,n))
   res%mode=opt%par; res%log_posterior=opt%value; res%iterations=opt%iterations; res%converged=opt%converged
   call numerical_hessian(f,res%mode,res%hessian); neg_h=-res%hessian
   call make_positive_definite(neg_h,info=info)
   call inverse_spd(neg_h,res%covariance,info)
   if(info/=0) then; res%covariance=0.0_dp; res%converged=.false.; end if
end subroutine laplace_approximation
end module ld_optim
