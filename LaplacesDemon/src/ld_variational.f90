module ld_variational
use ld_kinds, only: dp
use ld_interfaces, only: log_target_iface
use ld_random, only: rand_uniform, rand_mvn
use ld_linalg, only: inverse_spd, make_positive_definite
use ld_numerics, only: numerical_gradient, numerical_hessian
implicit none
private
public :: vb_result_t, variational_bayes_salimans2

type :: vb_result_t
   real(dp), allocatable :: mean(:), covariance(:,:), history_mean(:,:), history_var(:,:)
   real(dp) :: log_target = -huge(1.0_dp)
   real(dp) :: tolerance = huge(1.0_dp)
   real(dp) :: step_size = 0.0_dp
   integer :: iterations = 0
   logical :: converged = .false.
end type vb_result_t

contains

subroutine variational_bayes_salimans2(f,x0,res,initial_cov,iterations,interval,stop_tolerance)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(vb_result_t), intent(out) :: res
   real(dp), intent(in), optional :: initial_cov(:,:),interval,stop_tolerance
   integer, intent(in), optional :: iterations
   integer :: p,maxit,it,info,tries,used
   real(dp) :: hstep,tol,w,half2,lp_old,lp_new
   real(dp), allocatable :: m(:),v(:,:),prec(:,:),a(:),z(:),xstar(:),g(:),h(:,:)
   real(dp), allocatable :: abar(:),zbar(:),pbar(:,:),vbar(:,:),mbar(:),mbar_last(:),vbar_last(:,:)
   real(dp), allocatable :: trial(:),histm(:,:),histv(:,:)
   p=size(x0); maxit=1000; if(present(iterations)) maxit=max(10,iterations)
   hstep=1.0e-6_dp; if(present(interval)) hstep=interval
   tol=1.0e-5_dp; if(present(stop_tolerance)) tol=stop_tolerance
   allocate(m(p),v(p,p),prec(p,p),a(p),z(p),xstar(p),g(p),h(p,p))
   allocate(abar(p),zbar(p),pbar(p,p),vbar(p,p),mbar(p),mbar_last(p),vbar_last(p,p),trial(p))
   allocate(histm(maxit,p),histv(maxit,p))
   m=x0; v=identity_local(p)
   if(present(initial_cov)) v=initial_cov
   call make_positive_definite(v); call inverse_spd(v,prec,info)
   if(info/=0) then; v=identity_local(p); prec=v; end if
   a=0.0_dp; z=m; abar=0.0_dp; zbar=0.0_dp; pbar=0.0_dp
   vbar=v; mbar=m; mbar_last=m; vbar_last=v
   w=1.0_dp/sqrt(real(maxit,dp)); half2=2.0_dp/real(maxit,dp)
   lp_old=f(m); used=maxit; res%tolerance=huge(1.0_dp)
   do it=1,maxit
      xstar=m; tries=0
      do
         call rand_mvn(m,v,trial,info); tries=tries+1
         if(info==0 .and. finite_local(f(trial))) then; xstar=trial; exit; end if
         if(tries>=20) exit
      end do
      call numerical_gradient(f,xstar,g,hstep)
      call numerical_hessian(f,xstar,h,max(1.0e-4_dp,100.0_dp*hstep))
      a=(1.0_dp-w)*a+w*g
      prec=(1.0_dp-w)*prec-w*h
      z=(1.0_dp-w)*z+w*xstar
      call make_positive_definite(prec); call inverse_spd(prec,v,info)
      if(info/=0) then; prec=identity_local(p); v=prec; end if
      trial=matmul(v,a)+z; lp_new=f(trial)
      if(finite_local(lp_new) .and. log(rand_uniform())<lp_new-lp_old) then
         m=trial; lp_old=lp_new
      end if
      histm(it,:)=m; histv(it,:)=diag_local(v)
      if(it>maxit/2) then
         mbar_last=mbar; vbar_last=vbar
         abar=abar+half2*g; pbar=pbar-half2*h; zbar=zbar+half2*xstar
         call make_positive_definite(pbar); call inverse_spd(pbar,vbar,info)
         if(info/=0) then; pbar=prec; vbar=v; end if
         mbar=matmul(vbar,abar)+zbar
         res%tolerance=sqrt(sum((mbar-mbar_last)**2))+sqrt(sum((diag_local(vbar)-diag_local(vbar_last))**2))
         if(res%tolerance<=tol) then; used=it; exit; end if
      end if
   end do
   if(used<=maxit/2) then; mbar=m; vbar=v; end if
   if(any(abs(m-mbar)>3.0_dp*sqrt(max(diag_local(vbar),1.0e-14_dp)))) then
      mbar=m; vbar=diag_matrix(diag_local(v))
   end if
   allocate(res%mean(p),res%covariance(p,p),res%history_mean(used,p),res%history_var(used,p))
   res%mean=mbar; res%covariance=vbar; res%history_mean=histm(1:used,:); res%history_var=histv(1:used,:)
   res%log_target=f(mbar); res%iterations=used; res%step_size=w; res%converged=(res%tolerance<=tol)
end subroutine variational_bayes_salimans2

pure logical function finite_local(x)
   real(dp), intent(in) :: x
   finite_local=(x==x .and. abs(x)<huge(1.0_dp))
end function finite_local

pure function identity_local(n) result(a)
   integer, intent(in) :: n
   real(dp) :: a(n,n)
   integer :: i
   a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
end function identity_local

pure function diag_local(a) result(d)
   real(dp), intent(in) :: a(:,:)
   real(dp) :: d(min(size(a,1),size(a,2)))
   integer :: i
   do i=1,size(d); d(i)=a(i,i); end do
end function diag_local

pure function diag_matrix(d) result(a)
   real(dp), intent(in) :: d(:)
   real(dp) :: a(size(d),size(d))
   integer :: i
   a=0.0_dp; do i=1,size(d); a(i,i)=d(i); end do
end function diag_matrix

end module ld_variational
