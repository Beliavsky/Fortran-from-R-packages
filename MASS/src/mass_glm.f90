! SPDX-License-Identifier: GPL-3.0-only
module mass_glm
  use rrcov_kinds, only : dp
  use rrcov_linalg, only : general_inverse
  use mass_types, only : regression_result, loglinear_result, mass_success, &
    mass_invalid_argument, mass_dimension_error, mass_no_convergence
  use mass_math, only : weighted_least_squares
  use mass_distribution, only : theta_ml, negative_binomial_logpmf
  implicit none
  private
  public :: poisson_glm_fit, glm_nb_fit, loglinear_fit
  public :: negative_binomial_variance, negative_binomial_deviance
contains

  pure elemental function negative_binomial_variance(mu,theta) result(value)
    real(dp),intent(in)::mu,theta
    real(dp)::value
    value=mu+mu*mu/theta
  end function negative_binomial_variance

  pure elemental function negative_binomial_deviance(y,mu,theta) result(value)
    real(dp),intent(in)::y,mu,theta
    real(dp)::value,part1
    if(y<=0.0_dp)then
      part1=0.0_dp
    else
      part1=y*log(y/mu)
    end if
    value=2.0_dp*(part1-(y+theta)*log((y+theta)/(mu+theta)))
  end function negative_binomial_deviance

  subroutine poisson_glm_fit(x,y,result,case_weights,offset,maxit,tolerance)
    real(dp),intent(in)::x(:,:),y(:)
    type(regression_result),intent(out)::result
    real(dp),intent(in),optional::case_weights(:),offset(:),tolerance
    integer,intent(in),optional::maxit
    real(dp),allocatable::cw(:),off(:),beta(:),oldbeta(:),eta(:),mu(:),z(:),w(:),res(:),cov(:,:)
    real(dp)::tol,delta,ll
    integer::n,p,it,mit,rank,status,i
    n=size(y);p=size(x,2);mit=100;if(present(maxit))mit=maxit
    tol=1.0e-8_dp;if(present(tolerance))tol=tolerance
    if(size(x,1)/=n .or. n<p .or. any(y<0.0_dp))then;result%status=mass_dimension_error;return;end if
    allocate(cw(n),off(n));cw=1.0_dp;off=0.0_dp
    if(present(case_weights))then
      if(size(case_weights)/=n .or. any(case_weights<0.0_dp))then;result%status=mass_invalid_argument;return;end if
      cw=case_weights
    end if
    if(present(offset))then
      if(size(offset)/=n)then;result%status=mass_invalid_argument;return;end if
      off=offset
    end if
    allocate(beta(p),oldbeta(p),eta(n),mu(n),z(n),w(n))
    beta=0.0_dp
    if(p>=1)beta(1)=log(max(sum(cw*y)/max(sum(cw),tiny(1.0_dp)),0.1_dp))
    do it=1,mit
      oldbeta=beta;eta=off+matmul(x,beta);mu=exp(min(eta,700.0_dp));w=cw*mu
      z=eta+(y-mu)/max(mu,tiny(1.0_dp))-off
      call weighted_least_squares(x,z,w,beta,res,cov,rank,status)
      delta=maxval(abs(beta-oldbeta))/max(1.0_dp,maxval(abs(oldbeta)))
      if(delta<tol)exit
    end do
    eta=off+matmul(x,beta);mu=exp(min(eta,700.0_dp));res=y-mu;ll=0.0_dp
    do i=1,n;ll=ll+cw(i)*(y(i)*log(max(mu(i),tiny(1.0_dp)))-mu(i)-log_gamma(y(i)+1.0_dp));end do
    result%coefficients=beta;result%fitted=mu;result%residuals=res;result%weights=cw*mu
    result%covariance=cov;result%rank=rank;result%df_residual=n-rank;result%iterations=it
    result%log_likelihood=ll;result%aic=-2.0_dp*ll+2.0_dp*real(rank,dp)
    result%sigma=1.0_dp;result%status=merge(mass_success,mass_no_convergence,it<=mit);result%method="Poisson GLM"
    allocate(result%leverages(n));result%leverages=0.0_dp
  end subroutine poisson_glm_fit

  subroutine glm_nb_fit(x,y,result,init_theta,case_weights,offset,maxit,tolerance)
    real(dp),intent(in)::x(:,:),y(:)
    type(regression_result),intent(out)::result
    real(dp),intent(in),optional::init_theta,case_weights(:),offset(:),tolerance
    integer,intent(in),optional::maxit
    real(dp),allocatable::cw(:),off(:),beta(:),oldbeta(:),eta(:),mu(:),z(:),w(:),res(:),cov(:,:)
    real(dp)::theta,se_theta,tol,delta,ll
    integer::n,p,it,mit,rank,status,i,theta_status
    n=size(y);p=size(x,2);mit=50;if(present(maxit))mit=maxit
    tol=1.0e-7_dp;if(present(tolerance))tol=tolerance
    if(size(x,1)/=n .or. n<=p .or. any(y<0.0_dp))then;result%status=mass_dimension_error;return;end if
    allocate(cw(n),off(n));cw=1.0_dp;off=0.0_dp
    if(present(case_weights))then
      if(size(case_weights)/=n.or.any(case_weights<0.0_dp))then;result%status=mass_invalid_argument;return;end if
      cw=case_weights
    end if
    if(present(offset))then
      if(size(offset)/=n)then;result%status=mass_invalid_argument;return;end if
      off=offset
    end if
    theta=10.0_dp;if(present(init_theta))theta=max(init_theta,1.0e-6_dp)
    allocate(beta(p),oldbeta(p),eta(n),mu(n),z(n),w(n));beta=0.0_dp
    if(p>=1)beta(1)=log(max(sum(cw*y)/max(sum(cw),tiny(1.0_dp)),0.1_dp))
    do it=1,mit
      oldbeta=beta
      eta=off+matmul(x,beta);mu=exp(min(eta,700.0_dp))
      w=cw*mu/(1.0_dp+mu/theta)
      z=eta+(y-mu)/max(mu,tiny(1.0_dp))-off
      call weighted_least_squares(x,z,w,beta,res,cov,rank,status)
      eta=off+matmul(x,beta);mu=exp(min(eta,700.0_dp))
      call theta_ml(y,mu,theta,se_theta,theta_status,weights=cw,limit=20,tolerance=tol)
      delta=maxval(abs(beta-oldbeta))/max(1.0_dp,maxval(abs(oldbeta)))
      if(delta<tol)exit
    end do
    res=y-mu;ll=0.0_dp
    do i=1,n;ll=ll+cw(i)*negative_binomial_logpmf(y(i),mu(i),theta);end do
    result%coefficients=beta;result%fitted=mu;result%residuals=res;result%weights=w
    result%covariance=cov;result%rank=rank;result%df_residual=n-rank;result%iterations=it
    result%log_likelihood=ll;result%aic=-2.0_dp*ll+2.0_dp*real(rank+1,dp);result%theta=theta
    result%sigma=1.0_dp;result%status=merge(mass_success,mass_no_convergence,it<=mit);result%method="negative binomial GLM"
    allocate(result%leverages(n));result%leverages=0.0_dp
  end subroutine glm_nb_fit

  subroutine loglinear_fit(design,counts,result,weights,maxit,tolerance)
    real(dp),intent(in)::design(:,:),counts(:)
    type(loglinear_result),intent(out)::result
    real(dp),intent(in),optional::weights(:),tolerance
    integer,intent(in),optional::maxit
    type(regression_result)::fit
    real(dp),allocatable::cw(:)
    real(dp)::term
    integer::i,n
    n=size(counts);allocate(cw(n));cw=1.0_dp;if(present(weights))cw=weights
    call poisson_glm_fit(design,counts,fit,cw,maxit=maxit,tolerance=tolerance)
    result%fitted=fit%fitted;result%coefficients=fit%coefficients;result%iterations=fit%iterations
    result%df_residual=fit%df_residual;result%status=fit%status;result%deviance=0.0_dp;result%pearson=0.0_dp
    do i=1,n
      if(counts(i)>0.0_dp)then;term=counts(i)*log(counts(i)/max(fit%fitted(i),tiny(1.0_dp)))
      else;term=0.0_dp;end if
      result%deviance=result%deviance+2.0_dp*cw(i)*(term-counts(i)+fit%fitted(i))
      result%pearson=result%pearson+cw(i)*(counts(i)-fit%fitted(i))**2/max(fit%fitted(i),tiny(1.0_dp))
    end do
  end subroutine loglinear_fit

end module mass_glm
