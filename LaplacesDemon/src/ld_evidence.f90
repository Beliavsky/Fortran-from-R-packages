module ld_evidence
use ld_kinds, only: dp, pi
use ld_interfaces, only: log_target_iface
use ld_linalg, only: chol_lower, weighted_mean, weighted_covariance, make_positive_definite, logdet_spd
use ld_quadrature, only: gauss_hermite_rule, gauss_hermite_cube
use ld_numerics, only: logsumexp
implicit none
private
public :: iterative_quad_result_t, iterative_gauss_hermite
public :: lml_laplace, lml_harmonic, lml_generalized_harmonic

type :: iterative_quad_result_t
   real(dp), allocatable :: mean(:), covariance(:,:)
   real(dp) :: log_normalizer=-huge(1.0_dp)
   integer :: iterations=0
   logical :: converged=.false.
end type
contains
function lml_laplace(logpost_at_mode,covariance) result(v)
   real(dp),intent(in)::logpost_at_mode,covariance(:,:)
   real(dp)::v
   v=0.5_dp*real(size(covariance,1),dp)*log(2.0_dp*pi)+0.5_dp*logdet_spd(covariance)+logpost_at_mode
end function lml_laplace

function lml_harmonic(log_lik) result(v)
   real(dp),intent(in)::log_lik(:); real(dp)::v,m
   m=maxval(-log_lik); v=-(m+log(sum(exp(-log_lik-m))/real(size(log_lik),dp)))
end function lml_harmonic

function lml_generalized_harmonic(log_lik,log_g) result(v)
   real(dp),intent(in)::log_lik(:),log_g(:); real(dp)::v,m
   real(dp)::x(size(log_lik))
   x=log_g-log_lik; m=maxval(x); v=-(m+log(sum(exp(x-m))/real(size(x),dp)))
end function lml_generalized_harmonic

subroutine iterative_gauss_hermite(f,mean0,cov0,order,max_iter,tol,res)
   procedure(log_target_iface)::f
   real(dp),intent(in)::mean0(:),cov0(:,:)
   integer,intent(in)::order,max_iter
   real(dp),intent(in)::tol
   type(iterative_quad_result_t),intent(out)::res
   integer::p,it,npts,i,info
   real(dp),allocatable::nodes1(:),w1(:),znodes(:,:),zw(:),theta(:,:),lw(:),w(:),l(:,:),newm(:),newc(:,:)
   real(dp)::lse,change,log_jac
   p=size(mean0)
   if(p>5) error stop 'iterative_gauss_hermite: dimensions > 5 are intentionally unsupported by the tensor rule'
   allocate(nodes1(order),w1(order)); call gauss_hermite_rule(order,nodes1,w1); call gauss_hermite_cube(nodes1,w1,p,znodes,zw)
   npts=size(zw); allocate(theta(npts,p),lw(npts),w(npts),l(p,p),newm(p),newc(p,p),res%mean(p),res%covariance(p,p))
   res%mean=mean0; res%covariance=cov0
   do it=1,max_iter
      call chol_lower(res%covariance,l,info); if(info/=0) exit
      do i=1,npts
         theta(i,:)=res%mean+sqrt(2.0_dp)*matmul(l,znodes(i,:))
         ! Hermite rule integrates exp(-z^2); undo that weight in the integrand.
         lw(i)=log(max(zw(i),tiny(1.0_dp)))+f(theta(i,:))+sum(znodes(i,:)**2)
      end do
      lse=logsumexp(lw); w=exp(lw-lse); newm=weighted_mean(theta,w); newc=weighted_covariance(theta,w,newm)
      newc=newc+1e-10_dp*identity(p); call make_positive_definite(newc)
      change=max(maxval(abs(newm-res%mean)),maxval(abs(newc-res%covariance)))
      log_jac=0.5_dp*real(p,dp)*log(2.0_dp)+sum(log([(l(i,i),i=1,p)]))
      res%mean=newm; res%covariance=newc; res%iterations=it; res%log_normalizer=lse+log_jac
      if(change<tol) then; res%converged=.true.; exit; end if
   end do
contains
   pure function identity(n) result(a)
      integer,intent(in)::n; real(dp)::a(n,n); integer::j
      a=0.0_dp; do j=1,n; a(j,j)=1.0_dp; end do
   end function identity
end subroutine iterative_gauss_hermite
end module ld_evidence
