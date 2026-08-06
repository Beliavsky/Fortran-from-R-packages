! SPDX-License-Identifier: GPL-3.0-only
module mass_robust
  use rrcov_kinds, only : dp
  use rrcov_types, only : covariance_result, rrcov_success
  use rrcov_robust, only : cov_mcd_rr => cov_mcd, cov_mve_rr => cov_mve, cov_classic, robust_covariance
  use rrcov_linalg, only : general_inverse, mahalanobis_squared
  use rrcov_random, only : seed_random, random_subset
  use mass_types, only : regression_result, mass_success, mass_invalid_argument, mass_dimension_error, mass_no_convergence
  use mass_math, only : least_squares, type7_quantile, mad_mass, sort_real_mass
  implicit none
  private
  public :: cov_trob, cov_mcd, cov_mve, cov_rob
  public :: lqs_fit, lmsreg, ltsreg
contains

  subroutine cov_trob(x,result,weights,center,nu,maxit,tolerance)
    real(dp),intent(in)::x(:,:)
    type(covariance_result),intent(out)::result
    real(dp),intent(in),optional::weights(:),center(:),nu,tolerance
    integer,intent(in),optional::maxit
    real(dp),allocatable::wt(:),w(:),w0(:),loc(:),xc(:,:),cov(:,:),inv(:,:),q(:)
    real(dp)::df,tol,sumwt
    integer::n,p,it,mit,i,st
    n=size(x,1);p=size(x,2);df=5.0_dp;if(present(nu))df=nu;tol=0.01_dp;if(present(tolerance))tol=tolerance
    mit=25;if(present(maxit))mit=maxit
    if(n<=p .or. p<1 .or. df<=0.0_dp)then;result%status=mass_invalid_argument;return;end if
    allocate(wt(n));wt=1.0_dp;if(present(weights))then
      if(size(weights)/=n .or. any(weights<0.0_dp) .or. sum(weights)<=0.0_dp)then;result%status=mass_invalid_argument;return;end if
      wt=weights
    end if
    sumwt=sum(wt);allocate(loc(p));loc=matmul(wt,x)/sumwt;if(present(center))then
      if(size(center)/=p)then;result%status=mass_dimension_error;return;end if
      loc=center
    end if
    allocate(w(n),w0(n),xc(n,p),cov(p,p),q(n));w=wt*(1.0_dp+real(p,dp)/df)
    do it=1,mit
      w0=w;xc=x-spread(loc,1,n);cov=0.0_dp
      do i=1,n;cov=cov+w(i)*outer(xc(i,:),xc(i,:));end do
      cov=cov/max(sum(w),tiny(1.0_dp));inv=general_inverse(cov,st)
      do i=1,n;q(i)=dot_product(xc(i,:),matmul(inv,xc(i,:)));end do
      w=wt*(df+real(p,dp))/(df+q)
      if(.not.present(center))loc=matmul(w,x)/sum(w)
      if(maxval(abs(w-w0))<tol)exit
    end do
    xc=x-spread(loc,1,n);cov=0.0_dp
    do i=1,n;cov=cov+w(i)*outer(xc(i,:),xc(i,:));end do
    cov=cov/sumwt
    result%center=loc;result%covariance=0.5_dp*(cov+transpose(cov));result%weights=w;result%distances=q
    allocate(result%subset(n));result%subset=[(i,i=1,n)];result%n_obs=n;result%rank=p;result%iterations=it
    result%status=merge(rrcov_success,4,it<=mit);result%method="Student-t robust covariance"
  contains
    pure function outer(a,b) result(c)
      real(dp),intent(in)::a(:),b(:);real(dp)::c(size(a),size(b));integer::j
      do j=1,size(a);c(j,:)=a(j)*b;end do
    end function outer
  end subroutine cov_trob

  subroutine cov_mcd(x,result,alpha,nsamp,seed,reweight)
    real(dp),intent(in)::x(:,:)
    type(covariance_result),intent(out)::result
    real(dp),intent(in),optional::alpha
    integer,intent(in),optional::nsamp,seed
    logical,intent(in),optional::reweight
    call cov_mcd_rr(x,result,alpha,nsamp,seed,reweight=reweight)
  end subroutine cov_mcd

  subroutine cov_mve(x,result,alpha,nsamp,seed,reweight)
    real(dp),intent(in)::x(:,:)
    type(covariance_result),intent(out)::result
    real(dp),intent(in),optional::alpha
    integer,intent(in),optional::nsamp,seed
    logical,intent(in),optional::reweight
    call cov_mve_rr(x,result,alpha,nsamp,seed,reweight)
  end subroutine cov_mve

  subroutine cov_rob(x,result,method,alpha,nsamp,seed)
    real(dp),intent(in)::x(:,:)
    type(covariance_result),intent(out)::result
    character(len=*),intent(in),optional::method
    real(dp),intent(in),optional::alpha
    integer,intent(in),optional::nsamp,seed
    character(len=16)::m
    m="mve";if(present(method))m=method
    call robust_covariance(x,m,result,alpha=alpha,nsamp=nsamp,seed=seed)
  end subroutine cov_rob

  subroutine lqs_fit(x,y,result,method,quantile_used,nsamp,seed,intercept)
    real(dp),intent(in)::x(:,:),y(:)
    type(regression_result),intent(out)::result
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::quantile_used,nsamp,seed
    logical,intent(in),optional::intercept
    real(dp),allocatable::design(:,:),beta(:),res(:),subset_res(:),refit_res(:)
    real(dp),allocatable::bestbeta(:),squared(:),work(:),cov(:,:)
    integer,allocatable::idx(:),inlier_idx(:)
    character(len=8)::m
    logical::add_intercept
    integer::n,p,q,trials,it,rank,status,st,i,h
    real(dp)::criterion,bestcrit,scale
    n=size(y);add_intercept=.true.;if(present(intercept))add_intercept=intercept
    if(size(x,1)/=n .or. n<2)then;result%status=mass_dimension_error;return;end if
    if(add_intercept)then;allocate(design(n,size(x,2)+1));design(:,1)=1.0_dp;design(:,2:)=x
    else;design=x;end if
    p=size(design,2);if(n<=p)then;result%status=mass_invalid_argument;return;end if
    m="lts";if(present(method))m=method
    q=(n+p+1)/2;if(present(quantile_used))q=max(p,min(n,quantile_used))
    trials=min(500,max(100,20*p));if(present(nsamp))trials=max(1,nsamp)
    call seed_random(seed)
    allocate(idx(p),bestbeta(p),squared(n),work(n),res(n))
    bestcrit=huge(1.0_dp)
    do it=1,trials
      if(it==1)then;idx=[(i,i=1,p)];else;call random_subset(n,p,idx);end if
      call least_squares(design(idx,:),y(idx),beta,subset_res,rank,status)
      if(status/=mass_success)cycle
      res=y-matmul(design,beta);squared=res*res;work=squared;call sort_real_mass(work)
      select case(trim(m))
      case("lms")
        criterion=type7_quantile(squared,0.5_dp)
      case default
        criterion=sum(work(1:q))
      end select
      if(criterion<bestcrit)then;bestcrit=criterion;bestbeta=beta;end if
    end do
    res=y-matmul(design,bestbeta);squared=res*res
    if(trim(m)=="lms")then;scale=1.4826_dp*(1.0_dp+5.0_dp/real(max(1,n-p),dp))*sqrt(type7_quantile(squared,0.5_dp))
    else;work=squared;call sort_real_mass(work);scale=sqrt(sum(work(1:q))/real(q,dp));end if
    h=count(abs(res)<=2.5_dp*max(scale,tiny(1.0_dp)))
    if(h>p)then
      inlier_idx=pack([(i,i=1,n)],abs(res)<=2.5_dp*max(scale,tiny(1.0_dp)))
      call least_squares(design(inlier_idx,:),y(inlier_idx),beta,refit_res,rank,status)
      bestbeta=beta;res=y-matmul(design,bestbeta)
    end if
    result%coefficients=bestbeta;result%residuals=res;result%fitted=y-res;result%sigma=scale
    result%rank=p;result%df_residual=n-p;result%iterations=trials;result%status=mass_success
    allocate(result%weights(n),result%leverages(n));result%weights=merge(1.0_dp,0.0_dp,abs(res)<=2.5_dp*max(scale,tiny(1.0_dp)))
    result%leverages=0.0_dp;result%method="least quantile squares "//trim(m)
    cov=general_inverse(matmul(transpose(design),design),st);result%covariance=scale*scale*cov
  end subroutine lqs_fit

  subroutine lmsreg(x,y,result,nsamp,seed,intercept)
    real(dp),intent(in)::x(:,:),y(:);type(regression_result),intent(out)::result
    integer,intent(in),optional::nsamp,seed;logical,intent(in),optional::intercept
    call lqs_fit(x,y,result,"lms",nsamp=nsamp,seed=seed,intercept=intercept)
  end subroutine lmsreg

  subroutine ltsreg(x,y,result,quantile_used,nsamp,seed,intercept)
    real(dp),intent(in)::x(:,:),y(:);type(regression_result),intent(out)::result
    integer,intent(in),optional::quantile_used,nsamp,seed;logical,intent(in),optional::intercept
    call lqs_fit(x,y,result,"lts",quantile_used,nsamp,seed,intercept)
  end subroutine ltsreg

end module mass_robust
