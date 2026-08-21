! SPDX-License-Identifier: GPL-2.0-or-later
module mnb_core
  use mnb_kinds, only : dp
  use mnb_types, only : mnb_fit_result
  use mnb_math, only : numerical_hessian, covariance_from_hessian, normal_cdf
  use mnb_optimizer, only : bfgs_maximize
  implicit none
  private
  public :: mnb_loglik, fit_mnb, mnb_cluster_sums
contains

  subroutine mnb_cluster_sums(y,mu,n,mi,ysum,musum)
    real(dp),intent(in)::y(:),mu(:)
    integer,intent(in)::n,mi
    real(dp),intent(out)::ysum(n),musum(n)
    integer::i,l,u
    do i=1,n
      l=(i-1)*mi+1;u=i*mi
      ysum(i)=sum(y(l:u));musum(i)=sum(mu(l:u))
    end do
  end subroutine mnb_cluster_sums

  real(dp) function mnb_loglik(par,y,x,n,mi,offset) result(ll)
    real(dp),intent(in)::par(:),y(:),x(:,:)
    integer,intent(in)::n,mi
    real(dp),intent(in),optional::offset(:)
    real(dp)::phi
    real(dp),allocatable::eta(:),mu(:),ys(:),ms(:)
    integer::i
    if(size(y)/=n*mi .or. size(x,1)/=size(y) .or. size(par)/=size(x,2)+1)then
      ll=-huge(1.0_dp);return
    end if
    phi=par(1)
    if(phi<=0.0_dp)then;ll=-huge(1.0_dp)/8.0_dp;return;end if
    allocate(eta(size(y)),mu(size(y)),ys(n),ms(n));eta=matmul(x,par(2:))
    if(present(offset))eta=eta+offset
    if(maxval(eta)>700.0_dp)then;ll=-huge(1.0_dp)/8.0_dp;return;end if
    mu=exp(eta);call mnb_cluster_sums(y,mu,n,mi,ys,ms)
    ll=0.0_dp
    do i=1,n
      ll=ll+log_gamma(phi+ys(i))-log_gamma(phi)+phi*log(phi)-(phi+ys(i))*log(phi+ms(i))
    end do
    ll=ll-sum(log_gamma(y+1.0_dp))+dot_product(y,eta)
  end function mnb_loglik

  function fit_mnb(start,y,x,n,mi,offset,maxit,tol) result(res)
    real(dp),intent(in)::start(:),y(:),x(:,:)
    integer,intent(in)::n,mi
    real(dp),intent(in),optional::offset(:)
    integer,intent(in),optional::maxit
    real(dp),intent(in),optional::tol
    type(mnb_fit_result)::res
    integer::mx,i
    real(dp)::tt
    logical::ok
    mx=500;if(present(maxit))mx=maxit;tt=1.0e-6_dp;if(present(tol))tt=tol
    allocate(res%par(size(start)),res%hessian(size(start),size(start)))
    allocate(res%covariance(size(start),size(start)),res%se(size(start)),res%z(size(start)),res%p_value(size(start)))
    call bfgs_maximize(obj,start,res%par,res%loglik,res%iterations,res%convergence,mx,tt)
    call numerical_hessian(obj,res%par,res%hessian)
    call covariance_from_hessian(res%hessian,res%covariance,ok)
    if(ok)then
      do i=1,size(start)
        res%se(i)=sqrt(max(0.0_dp,res%covariance(i,i)))
        if(res%se(i)>0.0_dp)then
          res%z(i)=res%par(i)/res%se(i);res%p_value(i)=2.0_dp*(1.0_dp-normal_cdf(abs(res%z(i))))
        else
          res%z(i)=0.0_dp;res%p_value(i)=1.0_dp
        end if
      end do
    else
      res%covariance=0.0_dp;res%se=huge(1.0_dp);res%z=0.0_dp;res%p_value=1.0_dp
    end if
  contains
    real(dp) function obj(p) result(v)
      real(dp),intent(in)::p(:)
      if(present(offset))then
        v=mnb_loglik(p,y,x,n,mi,offset)
      else
        v=mnb_loglik(p,y,x,n,mi)
      end if
    end function
  end function fit_mnb
end module mnb_core
