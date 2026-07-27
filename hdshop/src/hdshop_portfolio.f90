! SPDX-License-Identifier: GPL-3.0-only
! Derived from HDShOP 0.1.7.
module hdshop_portfolio
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use hdshop_kinds, only: dp
  use hdshop_linalg, only: inverse_matrix, pseudo_inverse_symmetric, &
                           quadratic_form, symmetrize
  use hdshop_stats, only: row_means, sample_covariance, normal_quantile, &
                          chi_square1_survival
  implicit none
  private

  type, public :: portfolio_result
    real(dp), allocatable :: covariance(:,:), inverse_covariance(:,:)
    real(dp), allocatable :: means(:), sample_weights(:), weights(:)
    real(dp), allocatable :: weight_intervals(:,:)
    real(dp) :: alpha = 1.0_dp
    real(dp) :: variance = 0.0_dp
    real(dp) :: mean_return = 0.0_dp
    real(dp) :: sharpe = 0.0_dp
    logical :: ok = .false.
    character(len=200) :: message = ''
  end type portfolio_result

  type, public :: frontier_result
    real(dp), allocatable :: frontier_sd(:), frontier_return(:)
    real(dp), allocatable :: portfolio_sd(:), portfolio_return(:)
    logical :: ok = .false.
    character(len=200) :: message = ''
  end type frontier_result

  public :: q_matrix, s_value, r_gmv, v_gmv, v_portfolio
  public :: alpha_mv_lt, alpha_mv_gt, alpha_gmv_lt, alpha_gmv_gt
  public :: mean_variance_portfolio, traditional_portfolio
  public :: shrinkage_mv_portfolio, shrinkage_gmv_portfolio
  public :: mv_shrink_portfolio, bayesian_frontier

contains

  function q_matrix(inv_cov) result(q)
    real(dp), intent(in) :: inv_cov(:,:)
    real(dp), allocatable :: q(:,:), ones(:), v(:)
    real(dp) :: d
    integer :: p
    p=size(inv_cov,1);allocate(q(p,p),ones(p));ones=1.0_dp
    v=matmul(inv_cov,ones);d=dot_product(ones,v)
    q=inv_cov-spread(v,2,p)*spread(v,1,p)/d
    call symmetrize(q)
  end function q_matrix

  real(dp) function s_value(mu,covariance,use_pinv) result(value)
    real(dp),intent(in)::mu(:),covariance(:,:)
    logical,intent(in),optional::use_pinv
    real(dp),allocatable::inv(:,:),q(:,:)
    logical::ok,pinv
    pinv=.false.;if(present(use_pinv))pinv=use_pinv
    if(pinv)then;call pseudo_inverse_symmetric(covariance,inv,ok);else;call inverse_matrix(covariance,inv,ok);endif
    if(.not.ok)then;value=huge(1.0_dp);return;endif
    q=q_matrix(inv);value=quadratic_form(mu,q)
  end function s_value

  real(dp) function r_gmv(mu,covariance,use_pinv) result(value)
    real(dp),intent(in)::mu(:),covariance(:,:)
    logical,intent(in),optional::use_pinv
    real(dp),allocatable::inv(:,:),ones(:),v(:)
    logical::ok,pinv
    pinv=.false.;if(present(use_pinv))pinv=use_pinv
    if(pinv)then;call pseudo_inverse_symmetric(covariance,inv,ok);else;call inverse_matrix(covariance,inv,ok);endif
    if(.not.ok)then;value=huge(1.0_dp);return;endif
    allocate(ones(size(mu)));ones=1.0_dp;v=matmul(inv,ones)
    value=dot_product(v,mu)/sum(v)
  end function r_gmv

  real(dp) function v_gmv(covariance,use_pinv) result(value)
    real(dp),intent(in)::covariance(:,:)
    logical,intent(in),optional::use_pinv
    real(dp),allocatable::inv(:,:),ones(:)
    logical::ok,pinv
    pinv=.false.;if(present(use_pinv))pinv=use_pinv
    if(pinv)then;call pseudo_inverse_symmetric(covariance,inv,ok);else;call inverse_matrix(covariance,inv,ok);endif
    if(.not.ok)then;value=huge(1.0_dp);return;endif
    allocate(ones(size(covariance,1)));ones=1.0_dp
    value=1.0_dp/dot_product(ones,matmul(inv,ones))
  end function v_gmv

  pure real(dp) function v_portfolio(covariance,w) result(value)
    real(dp),intent(in)::covariance(:,:),w(:)
    value=quadratic_form(w,covariance)
  end function v_portfolio

  pure real(dp) function inv_gamma(gamma) result(value)
    real(dp),intent(in)::gamma
    if(ieee_is_finite(gamma))then
      value=1.0_dp/gamma
    else
      value=0.0_dp
    end if
  end function inv_gamma

  pure real(dp) function alpha_mv_lt(gamma,c,s,r_g,r_b,v_c,v_b) result(alpha)
    real(dp),intent(in)::gamma,c,s,r_g,r_b,v_c,v_b
    real(dp)::ig,num,den
    ig=inv_gamma(gamma)
    num=(r_g-r_b)*(1.0_dp+1.0_dp/(1.0_dp-c))*ig+(v_b-v_c)+s*ig*ig/(1.0_dp-c)
    den=v_c/(1.0_dp-c)-2.0_dp*(v_c+(r_b-r_g)*ig/(1.0_dp-c))+ &
      (s+c)*ig*ig/(1.0_dp-c)**3+v_b
    alpha=num/den
  end function alpha_mv_lt

  pure real(dp) function alpha_mv_gt(gamma,c,s,r_g,r_b,v_c,v_b) result(alpha)
    real(dp),intent(in)::gamma,c,s,r_g,r_b,v_c,v_b
    real(dp)::ig,num,den
    ig=inv_gamma(gamma)
    num=(r_g-r_b)*(1.0_dp+1.0_dp/(c*(c-1.0_dp)))*ig+(v_b-v_c)+ &
      s*ig*ig/(c*(c-1.0_dp))
    den=c*c*v_c/(c-1.0_dp)-2.0_dp*(v_c+(r_b-r_g)*ig/(c*(c-1.0_dp)))+ &
      (s+c*c)*ig*ig/(c-1.0_dp)**3+v_b
    alpha=num/den
  end function alpha_mv_gt

  pure real(dp) function alpha_gmv_lt(c,v_sample,v_b) result(alpha)
    real(dp),intent(in)::c,v_sample,v_b
    real(dp)::v_c,num
    v_c=v_sample/(1.0_dp-c);num=(1.0_dp-c)*(v_b-v_c)
    alpha=num/(num+c*v_c)
  end function alpha_gmv_lt

  pure real(dp) function alpha_gmv_gt(c,v_c,v_b) result(alpha)
    real(dp),intent(in)::c,v_c,v_b
    alpha=(v_b-v_c)/(c*c*v_c/(c-1.0_dp)-2.0_dp*v_c+v_b)
  end function alpha_gmv_gt

  function mean_variance_portfolio(mean_vec,covariance,gamma,use_pinv) result(res)
    real(dp),intent(in)::mean_vec(:),covariance(:,:),gamma
    logical,intent(in),optional::use_pinv
    type(portfolio_result)::res
    real(dp),allocatable::inv(:,:),ones(:),q(:,:),w(:)
    real(dp)::den,ig
    logical::ok,pinv
    integer::p
    p=size(mean_vec);pinv=.false.;if(present(use_pinv))pinv=use_pinv
    call allocate_result(res,p)
    if(any(shape(covariance)/=[p,p]))then;res%message='dimension mismatch';return;endif
    if(pinv)then;call pseudo_inverse_symmetric(covariance,inv,ok);else;call inverse_matrix(covariance,inv,ok);endif
    if(.not.ok)then;res%message='covariance matrix is singular';return;endif
    allocate(ones(p));ones=1.0_dp;q=q_matrix(inv);den=dot_product(ones,matmul(inv,ones));ig=inv_gamma(gamma)
    w=matmul(inv,ones)/den+ig*matmul(q,mean_vec)
    res%covariance=covariance;res%inverse_covariance=inv;res%means=mean_vec
    res%sample_weights=w;res%weights=w;call finish_result(res);res%ok=.true.
  end function mean_variance_portfolio

  function traditional_portfolio(x,gamma) result(res)
    real(dp),intent(in)::x(:,:),gamma
    type(portfolio_result)::res
    real(dp),allocatable::cov(:,:),mu(:)
    cov=sample_covariance(x);mu=row_means(x)
    res=mean_variance_portfolio(mu,cov,gamma,use_pinv=size(x,1)>=size(x,2))
  end function traditional_portfolio

  function shrinkage_mv_portfolio(x,gamma,target,beta) result(res)
    real(dp),intent(in)::x(:,:),gamma,target(:)
    real(dp),intent(in),optional::beta
    type(portfolio_result)::res
    real(dp),allocatable::cov(:,:),inv(:,:),mu(:),ones(:),q(:,:),w0(:),w(:)
    real(dp)::c,den,v_sample,v_c,v_b,r_g,r_b,s_c,al,level,z,eta,w_est,omega,td
    logical::ok
    integer::p,n,i
    p=size(x,1);n=size(x,2);call allocate_result(res,p)
    if(size(target)/=p .or. p==n)then;res%message='target mismatch or square data matrix';return;endif
    cov=sample_covariance(x);mu=row_means(x);allocate(ones(p));ones=1.0_dp;c=real(p,dp)/real(n,dp)
    if(p<n)then;call inverse_matrix(cov,inv,ok);else;call pseudo_inverse_symmetric(cov,inv,ok);endif
    if(.not.ok)then;res%message='covariance inversion failed';return;endif
    q=q_matrix(inv);den=dot_product(ones,matmul(inv,ones));w0=matmul(inv,ones)/den+inv_gamma(gamma)*matmul(q,mu)
    v_sample=1.0_dp/den;v_b=quadratic_form(target,cov);r_g=dot_product(matmul(inv,ones),mu)/den
    r_b=dot_product(target,mu)
    if(p<n)then
      v_c=v_sample/(1.0_dp-c);s_c=(1.0_dp-c)*quadratic_form(mu,q)-c
      al=alpha_mv_lt(gamma,c,s_c,r_g,r_b,v_c,v_b)
    else
      v_c=v_sample/(c*(c-1.0_dp));s_c=c*(c-1.0_dp)*quadratic_form(mu,q)-c
      al=alpha_mv_gt(gamma,c,s_c,r_g,r_b,v_c,v_b)
    end if
    w=al*w0+(1.0_dp-al)*target
    res%covariance=cov;res%inverse_covariance=inv;res%means=mu;res%sample_weights=w0;res%weights=w;res%alpha=al
    if(p<n)then
      allocate(res%weight_intervals(p,5));res%weight_intervals=0.0_dp
      level=0.05_dp;if(present(beta))level=beta;z=normal_quantile(1.0_dp-level/2.0_dp)
      do i=1,p
        if(abs(s_c)>sqrt(epsilon(1.0_dp)) .and. quadratic_form(mu,q)>0.0_dp)then
          eta=((s_c+c)/s_c)*dot_product(q(i,:),mu)/quadratic_form(mu,q)
        else
          eta=0.0_dp
        end if
        w_est=dot_product(inv(i,:),ones)/den+inv_gamma(gamma)*s_c*eta
        omega=(inv_gamma(gamma)**2*(s_c+1.0_dp)+v_c)*(1.0_dp-c)*q(i,i)+ &
          inv_gamma(gamma)**2*(s_c+c)**2*eta*eta
        omega=max(omega,0.0_dp)
        if(omega>0.0_dp)then
          td=real(n-p,dp)*w_est*w_est/omega
          res%weight_intervals(i,:)=[w(i),w_est-z*sqrt(omega/real(n-p,dp)), &
            w_est+z*sqrt(omega/real(n-p,dp)),td,chi_square1_survival(td)]
        else
          res%weight_intervals(i,:)=[w(i),w_est,w_est,0.0_dp,1.0_dp]
        end if
      end do
    end if
    call finish_result(res);res%ok=.true.
  end function shrinkage_mv_portfolio

  function shrinkage_gmv_portfolio(x,target,beta) result(res)
    real(dp),intent(in)::x(:,:),target(:)
    real(dp),intent(in),optional::beta
    type(portfolio_result)::res
    real(dp),allocatable::cov(:,:),inv(:,:),mu(:),ones(:),q(:,:),w0(:),w(:)
    real(dp)::c,den,v_sample,v_c,v_b,al,level,z,omega,td,w_est,lb
    logical::ok
    integer::p,n,i
    p=size(x,1);n=size(x,2);call allocate_result(res,p)
    if(size(target)/=p .or. p==n)then;res%message='target mismatch or square data matrix';return;endif
    cov=sample_covariance(x);mu=row_means(x);allocate(ones(p));ones=1.0_dp;c=real(p,dp)/real(n,dp)
    if(p<n)then;call inverse_matrix(cov,inv,ok);else;call pseudo_inverse_symmetric(cov,inv,ok);endif
    if(.not.ok)then;res%message='covariance inversion failed';return;endif
    den=dot_product(ones,matmul(inv,ones));w0=matmul(inv,ones)/den;v_sample=1.0_dp/den;v_b=quadratic_form(target,cov)
    if(p<n)then
      lb=(1.0_dp-c)*v_b*den-1.0_dp
      al=(1.0_dp-c)*lb/(c+(1.0_dp-c)*lb)
      v_c=v_sample/(1.0_dp-c)
    else
      v_c=v_sample/(c*(c-1.0_dp));al=alpha_gmv_gt(c,v_c,v_b)
    end if
    w=al*w0+(1.0_dp-al)*target
    res%covariance=cov;res%inverse_covariance=inv;res%means=mu;res%sample_weights=w0;res%weights=w;res%alpha=al
    if(p<n)then
      q=q_matrix(inv);allocate(res%weight_intervals(p,5));level=0.05_dp;if(present(beta))level=beta
      z=normal_quantile(1.0_dp-level/2.0_dp)
      do i=1,p
        w_est=w0(i);omega=max(v_c*(1.0_dp-c)*q(i,i),0.0_dp)
        if(omega>0.0_dp)then
          td=real(n-p,dp)*w_est*w_est/omega
          res%weight_intervals(i,:)=[w(i),w_est-z*sqrt(omega/real(n-p,dp)), &
            w_est+z*sqrt(omega/real(n-p,dp)),td,chi_square1_survival(td)]
        else
          res%weight_intervals(i,:)=[w(i),w_est,w_est,0.0_dp,1.0_dp]
        end if
      end do
    end if
    call finish_result(res);res%ok=.true.
  end function shrinkage_gmv_portfolio

  function mv_shrink_portfolio(x,gamma,kind,target,beta) result(res)
    real(dp),intent(in)::x(:,:),gamma
    character(len=*),intent(in)::kind
    real(dp),intent(in),optional::target(:),beta
    type(portfolio_result)::res
    if(trim(kind)=='traditional')then
      res=traditional_portfolio(x,gamma)
    else if(trim(kind)=='shrinkage')then
      if(.not.present(target))then
        call allocate_result(res,size(x,1));res%message='shrinkage requires target weights'
      else if(ieee_is_finite(gamma))then
        res=shrinkage_mv_portfolio(x,gamma,target,beta)
      else
        res=shrinkage_gmv_portfolio(x,target,beta)
      end if
    else
      call allocate_result(res,size(x,1));res%message='unknown portfolio kind'
    end if
  end function mv_shrink_portfolio

  function bayesian_frontier(x,portfolio_weights,npoints) result(res)
    real(dp),intent(in)::x(:,:),portfolio_weights(:,:)
    integer,intent(in),optional::npoints
    type(frontier_result)::res
    real(dp),allocatable::cov(:,:),mu(:),sigma_b(:,:),inv(:,:),ones(:),q(:,:)
    real(dp)::cb,vb,rb,curvature,maxsd,t
    logical::ok
    integer::p,n,m,k
    p=size(x,1);n=size(x,2);m=100;if(present(npoints))m=npoints
    allocate(res%frontier_sd(m),res%frontier_return(m))
    allocate(res%portfolio_sd(size(portfolio_weights,2)),res%portfolio_return(size(portfolio_weights,2)))
    if(size(portfolio_weights,1)/=p .or. n<=p+2)then;res%message='invalid weights or n <= p + 2';return;endif
    cov=sample_covariance(x);mu=row_means(x)
    do k=1,size(portfolio_weights,2)
      res%portfolio_sd(k)=sqrt(max(quadratic_form(portfolio_weights(:,k),cov),0.0_dp))
      res%portfolio_return(k)=dot_product(portfolio_weights(:,k),mu)
    end do
    sigma_b=real(n-1,dp)*cov;call inverse_matrix(sigma_b,inv,ok)
    if(.not.ok)then;res%message='singular Bayesian covariance';return;endif
    allocate(ones(p));ones=1.0_dp
    cb=1.0_dp/real(n-p-1,dp)+real(2*n-p-1,dp)/real(n*(n-p-1)*(n-p-2),dp)
    vb=cb/dot_product(ones,matmul(inv,ones));q=q_matrix(inv)
    rb=dot_product(ones,matmul(inv,mu))/dot_product(ones,matmul(inv,ones))
    curvature=quadratic_form(mu,q);maxsd=max(sqrt(vb),1.5_dp*maxval(res%portfolio_sd))
    do k=1,m
      if(m==1)then;t=0.0_dp;else;t=real(k-1,dp)/real(m-1,dp);endif
      res%frontier_sd(k)=sqrt(vb)+t*(maxsd-sqrt(vb))
      res%frontier_return(k)=rb+sqrt(max(curvature*(res%frontier_sd(k)**2-vb)/cb,0.0_dp))
    end do
    res%ok=.true.
  end function bayesian_frontier

  subroutine allocate_result(res,p)
    type(portfolio_result),intent(inout)::res
    integer,intent(in)::p
    allocate(res%covariance(p,p),res%inverse_covariance(p,p),res%means(p), &
      res%sample_weights(p),res%weights(p))
    res%covariance=0.0_dp;res%inverse_covariance=0.0_dp;res%means=0.0_dp
    res%sample_weights=0.0_dp;res%weights=0.0_dp
  end subroutine allocate_result

  subroutine finish_result(res)
    type(portfolio_result),intent(inout)::res
    res%variance=quadratic_form(res%weights,res%covariance)
    res%mean_return=dot_product(res%means,res%weights)
    if(res%variance>0.0_dp)res%sharpe=res%mean_return/sqrt(res%variance)
  end subroutine finish_result

end module hdshop_portfolio
