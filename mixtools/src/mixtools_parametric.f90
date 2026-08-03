! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_parametric
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mixtools_kinds, only : dp
  use mixtools_status
  use mixtools_types
  use mixtools_distributions
  use mixtools_linalg, only : inverse_spd, solve_least_squares
  use mixtools_utilities, only : sort_real
  implicit none
  private
  public :: normalmix_em, normalmix_em2comp, normalmix_mmlc, tauequivnormalmix_em
  public :: mvnormalmix_em, gammamix_em, multmix_em, repnormmix_em
contains
  subroutine initial_univariate(x,k,lambda,mu,sigma,status)
    real(dp),intent(in)::x(:)
    integer,intent(in)::k
    real(dp),allocatable,intent(out)::lambda(:),mu(:),sigma(:)
    integer,intent(out)::status
    real(dp),allocatable::xs(:)
    real(dp)::sd
    integer::j,n,idx
    n=size(x)
    if(n<k.or.k<1)then;allocate(lambda(0),mu(0),sigma(0));status=MIXTOOLS_INVALID_ARGUMENT;return;end if
    allocate(xs(n));xs=x;call sort_real(xs)
    allocate(lambda(k),mu(k),sigma(k));lambda=1.0_dp/real(k,dp)
    sd=sqrt(sum((x-sum(x)/real(n,dp))**2)/real(max(1,n-1),dp));sd=max(sd,1.0e-3_dp)
    do j=1,k
      idx=max(1,min(n,nint((real(j,dp)-0.5_dp)*real(n,dp)/real(k,dp))))
      mu(j)=xs(idx);sigma(j)=sd
    end do
    status=MIXTOOLS_SUCCESS
  end subroutine initial_univariate

  subroutine normal_estep(x,lambda,mu,sigma,post,ll,status)
    real(dp),intent(in)::x(:),lambda(:),mu(:),sigma(:)
    real(dp),intent(out)::post(size(x),size(lambda)),ll
    integer,intent(out)::status
    real(dp),allocatable::lw(:)
    real(dp)::ln
    integer::i,j,k
    k=size(lambda);allocate(lw(k));ll=0.0_dp
    do i=1,size(x)
      do j=1,k
        lw(j)=log(max(lambda(j),tiny(1.0_dp)))+normal_logpdf(x(i),mu(j),sigma(j))
      end do
      call normalize_logweights(lw,post(i,:),ln)
      if(.not.ieee_is_finite(ln))then;status=MIXTOOLS_NUMERICAL_ERROR;return;end if
      ll=ll+ln
    end do
    status=MIXTOOLS_SUCCESS
  end subroutine normal_estep

  subroutine normalmix_em(x,k,result,control,lambda0,mu0,sigma0,common_mean,common_sigma)
    real(dp),intent(in)::x(:)
    integer,intent(in)::k
    type(mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::lambda0(:),mu0(:),sigma0(:)
    logical,intent(in),optional::common_mean,common_sigma
    type(em_control)::ctl
    real(dp),allocatable::lambda(:),mu(:),sigma(:),post(:,:),history(:),nk(:)
    real(dp)::ll,newll,diff,global_mu,global_var
    logical::cm,cs
    integer::j,iter,status,n
    ctl=em_control();if(present(control))ctl=control
    cm=.false.;if(present(common_mean))cm=common_mean
    cs=.false.;if(present(common_sigma))cs=common_sigma
    n=size(x)
    call initial_univariate(x,k,lambda,mu,sigma,status)
    if(status/=MIXTOOLS_SUCCESS)then;result%status=status;return;end if
    if(present(lambda0))then;if(size(lambda0)==k)lambda=max(lambda0,tiny(1.0_dp));lambda=lambda/sum(lambda);end if
    if(present(mu0))then;if(size(mu0)==k)mu=mu0;end if
    if(present(sigma0))then;if(size(sigma0)==k)sigma=max(sigma0,ctl%minimum_scale);end if
    allocate(post(n,k),history(ctl%max_iterations+1),nk(k))
    call normal_estep(x,lambda,mu,sigma,post,ll,status)
    if(status/=MIXTOOLS_SUCCESS)then;result%status=status;return;end if
    history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1);lambda=max(nk/real(n,dp),tiny(1.0_dp));lambda=lambda/sum(lambda)
      if(cm)then
        global_mu=sum(x)/real(n,dp);mu=global_mu
      else
        do j=1,k;mu(j)=sum(post(:,j)*x)/max(nk(j),tiny(1.0_dp));end do
      end if
      if(cs)then
        global_var=0.0_dp
        do j=1,k;global_var=global_var+sum(post(:,j)*(x-mu(j))**2);end do
        sigma=sqrt(max(global_var/real(n,dp),ctl%minimum_scale**2))
      else
        do j=1,k
          sigma(j)=sqrt(max(sum(post(:,j)*(x-mu(j))**2)/max(nk(j),tiny(1.0_dp)),ctl%minimum_scale**2))
        end do
      end if
      call normal_estep(x,lambda,mu,sigma,post,newll,status)
      if(status/=MIXTOOLS_SUCCESS)exit
      history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=lambda;result%mu=mu;result%sigma=sigma;result%posterior=post
    result%iterations=min(iter,ctl%max_iterations);result%loglik=ll
    result%loglik_history=history(:result%iterations+1)
    result%converged=(iter<=ctl%max_iterations.and.status==MIXTOOLS_SUCCESS)
    result%status=merge(MIXTOOLS_SUCCESS,MIXTOOLS_NOT_CONVERGED,result%converged)
  end subroutine normalmix_em

  subroutine normalmix_em2comp(x,result,control,lambda0,mu0,variance0)
    real(dp),intent(in)::x(:)
    type(mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::lambda0,mu0(2),variance0(2)
    real(dp)::lam(2),mu(2),sig(2)
    lam=[0.5_dp,0.5_dp];if(present(lambda0))lam=[lambda0,1.0_dp-lambda0]
    mu=[minval(x),maxval(x)];if(present(mu0))mu=mu0
    sig=sqrt(max(sum((x-sum(x)/real(size(x),dp))**2)/real(max(1,size(x)-1),dp),1.0e-6_dp))
    if(present(variance0))sig=sqrt(max(variance0,1.0e-12_dp))
    call normalmix_em(x,2,result,control,lam,mu,sig)
  end subroutine normalmix_em2comp

  subroutine normalmix_mmlc(x,mean_matrix,mean_constant,var_matrix,gamma0,result,control)
    real(dp),intent(in)::x(:),mean_matrix(:,:),mean_constant(:),var_matrix(:,:),gamma0(:)
    type(mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    type(em_control)::ctl
    integer::k,p,q,n,j,iter,status
    real(dp),allocatable::lambda(:),mu(:),sigma(:),post(:,:),history(:),gamma(:),iv(:),nk(:)
    real(dp),allocatable::bmat(:,:),rhs(:),beta(:),cov(:,:),weights(:),r0(:,:),den(:),num(:)
    real(dp)::ll,newll,diff,rss
    ctl=em_control();if(present(control))ctl=control
    n=size(x);k=size(mean_matrix,1);p=size(mean_matrix,2);q=size(var_matrix,2)
    if(size(mean_constant)/=k.or.size(var_matrix,1)/=k.or.size(gamma0)/=q)then
      result%status=MIXTOOLS_DIMENSION_ERROR;return
    end if
    allocate(lambda(k),mu(k),sigma(k),post(n,k),history(ctl%max_iterations+1),gamma(q),iv(k),nk(k))
    lambda=1.0_dp/real(k,dp);gamma=max(gamma0,1.0e-6_dp);iv=matmul(var_matrix,gamma)
    if(any(iv<=0.0_dp))then;result%status=MIXTOOLS_INVALID_ARGUMENT;return;end if
    sigma=1.0_dp/sqrt(iv);mu=mean_constant
    call normal_estep(x,lambda,mu,sigma,post,ll,status);history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1);lambda=nk/real(n,dp)
      allocate(bmat(p,p),rhs(p),weights(k));bmat=0.0_dp;rhs=0.0_dp;weights=nk*iv
      do j=1,k
        bmat=bmat+weights(j)*outer(mean_matrix(j,:))
        rhs=rhs+mean_matrix(j,:)*iv(j)*sum(post(:,j)*(x-mean_constant(j)))
      end do
      call solve_least_squares(bmat,rhs,beta,cov,rss,status)
      if(status/=MIXTOOLS_SUCCESS)then
        call ridge_solve(bmat,rhs,beta,status)
        if(status/=MIXTOOLS_SUCCESS)exit
      end if
      mu=matmul(mean_matrix,beta)+mean_constant
      deallocate(bmat,rhs,weights)
      allocate(r0(n,k),den(q),num(q));do j=1,k;r0(:,j)=post(:,j)*(x-mu(j))**2;end do
      den=sum(matmul(r0,var_matrix),dim=1);num=0.0_dp
      do j=1,k;num=num+sum(post(:,j))*var_matrix(j,:)/iv(j);end do
      gamma=gamma*max(num,tiny(1.0_dp))/max(den,tiny(1.0_dp));iv=matmul(var_matrix,gamma)
      if(any(iv<=0.0_dp))then;status=MIXTOOLS_NUMERICAL_ERROR;exit;end if
      sigma=1.0_dp/sqrt(iv);deallocate(r0,den,num)
      call normal_estep(x,lambda,mu,sigma,post,newll,status);if(status/=MIXTOOLS_SUCCESS)exit
      history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=lambda;result%mu=mu;result%sigma=sigma;result%posterior=post;result%loglik=ll
    result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations.and.status==MIXTOOLS_SUCCESS
    result%status=merge(MIXTOOLS_SUCCESS,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    pure function outer(v) result(a)
      real(dp),intent(in)::v(:);real(dp)::a(size(v),size(v));integer::ii,jj
      do ii=1,size(v);do jj=1,size(v);a(ii,jj)=v(ii)*v(jj);end do;end do
    end function outer
    subroutine ridge_solve(a,b,x,stat)
      real(dp),intent(in)::a(:,:),b(:);real(dp),allocatable,intent(out)::x(:);integer,intent(out)::stat
      real(dp),allocatable::aa(:,:),c(:,:);real(dp)::r;integer::ii
      aa=a;do ii=1,size(a,1);aa(ii,ii)=aa(ii,ii)+1.0e-8_dp;end do
      call solve_least_squares(aa,b,x,c,r,stat)
    end subroutine ridge_solve
  end subroutine normalmix_mmlc

  subroutine tauequivnormalmix_em(x,result,control)
    real(dp),intent(in)::x(:)
    type(mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    real(dp)::m(3,2),c(3),a(3,2),g(2)
    m=reshape([1.0_dp,1.0_dp,1.0_dp,0.0_dp,-1.0_dp,1.0_dp],[3,2])
    c=0.0_dp;a=reshape([1.0_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp],[3,2]);g=[1.0_dp,1.0_dp]
    call normalmix_mmlc(x,m,c,a,g,result,control)
  end subroutine tauequivnormalmix_em

  subroutine mvnormalmix_em(x,k,result,control,common_mean,common_covariance)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::k
    type(mv_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    logical,intent(in),optional::common_mean,common_covariance
    type(em_control)::ctl
    integer::n,p,i,j,a,b,iter,status
    logical::cm,cc
    real(dp),allocatable::lambda(:),mu(:,:),sigma(:,:,:),post(:,:),history(:),lw(:),nk(:),d(:)
    real(dp)::ll,newll,ln,diff,global_mu(size(x,2)),pooled(size(x,2),size(x,2))
    ctl=em_control();if(present(control))ctl=control
    cm=.false.;if(present(common_mean))cm=common_mean
    cc=.false.;if(present(common_covariance))cc=common_covariance
    n=size(x,1);p=size(x,2)
    if(n<k.or.k<1)then;result%status=MIXTOOLS_INVALID_ARGUMENT;return;end if
    allocate(lambda(k),mu(p,k),sigma(p,p,k),post(n,k),history(ctl%max_iterations+1),lw(k),nk(k),d(p))
    lambda=1.0_dp/real(k,dp);global_mu=sum(x,dim=1)/real(n,dp);pooled=0.0_dp
    do i=1,n;d=x(i,:)-global_mu;do a=1,p;do b=1,p;pooled(a,b)=pooled(a,b)+d(a)*d(b);end do;end do;end do
    pooled=pooled/real(max(1,n-1),dp);do a=1,p;pooled(a,a)=pooled(a,a)+ctl%ridge;end do
    do j=1,k;mu(:,j)=x(1+mod((j-1)*max(1,n/k),n),:);sigma(:,:,j)=pooled;end do
    call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1);lambda=max(nk/real(n,dp),tiny(1.0_dp));lambda=lambda/sum(lambda)
      if(cm)then;mu=spread(sum(x,dim=1)/real(n,dp),2,k);else
        do j=1,k;mu(:,j)=matmul(transpose(x),post(:,j))/max(nk(j),tiny(1.0_dp));end do
      end if
      if(cc)then
        pooled=0.0_dp
        do j=1,k;do i=1,n;d=x(i,:)-mu(:,j);do a=1,p;do b=1,p
          pooled(a,b)=pooled(a,b)+post(i,j)*d(a)*d(b)
        end do;end do;end do;end do
        pooled=pooled/real(n,dp);do a=1,p;pooled(a,a)=pooled(a,a)+ctl%ridge;end do
        do j=1,k;sigma(:,:,j)=pooled;end do
      else
        do j=1,k
          sigma(:,:,j)=0.0_dp
          do i=1,n;d=x(i,:)-mu(:,j);do a=1,p;do b=1,p
            sigma(a,b,j)=sigma(a,b,j)+post(i,j)*d(a)*d(b)
          end do;end do;end do
          sigma(:,:,j)=sigma(:,:,j)/max(nk(j),tiny(1.0_dp));do a=1,p;sigma(a,a,j)=sigma(a,a,j)+ctl%ridge;end do
        end do
      end if
      call estep(newll);if(status/=MIXTOOLS_SUCCESS)exit
      history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=lambda;result%mu=mu;result%sigma=sigma;result%posterior=post;result%loglik=ll
    result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations.and.status==MIXTOOLS_SUCCESS
    result%status=merge(MIXTOOLS_SUCCESS,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      ll=0.0_dp;status=MIXTOOLS_SUCCESS
      do i=1,n
        do j=1,k;lw(j)=log(max(lambda(j),tiny(1.0_dp)))+logdmvnorm(x(i,:),mu(:,j),sigma(:,:,j),status);if(status/=0)return;end do
        call normalize_logweights(lw,post(i,:),ln);ll=ll+ln
      end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine mvnormalmix_em

  subroutine gammamix_em(x,k,result,control,fixed_shape,common_shape)
    real(dp),intent(in)::x(:)
    integer,intent(in)::k
    type(gamma_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::fixed_shape(:)
    logical,intent(in),optional::common_shape
    type(em_control)::ctl
    real(dp),allocatable::lambda(:),shape(:),scale(:),post(:,:),history(:),lw(:),nk(:)
    real(dp)::ll,newll,ln,diff,m,v,meanlog,target,a,newa,common
    integer::i,j,iter,status,n,step
    logical::cs
    ctl=em_control();if(present(control))ctl=control
    cs=.false.;if(present(common_shape))cs=common_shape;n=size(x)
    if(any(x<=0.0_dp).or.n<k)then;result%status=MIXTOOLS_INVALID_ARGUMENT;return;end if
    allocate(lambda(k),shape(k),scale(k),post(n,k),history(ctl%max_iterations+1),lw(k),nk(k))
    lambda=1.0_dp/real(k,dp);m=sum(x)/real(n,dp);v=sum((x-m)**2)/real(max(1,n-1),dp)
    shape=max(m*m/max(v,1.0e-6_dp),0.2_dp);scale=max(v/max(m,1.0e-6_dp),1.0e-3_dp)
    do j=1,k;scale(j)=scale(j)*(0.5_dp+real(j,dp)/real(k,dp));end do
    if(present(fixed_shape))then;if(size(fixed_shape)==1)shape=fixed_shape(1);if(size(fixed_shape)==k)shape=fixed_shape;end if
    call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1);lambda=nk/real(n,dp)
      do j=1,k
        m=sum(post(:,j)*x)/max(nk(j),tiny(1.0_dp));meanlog=sum(post(:,j)*log(x))/max(nk(j),tiny(1.0_dp))
        if(.not.present(fixed_shape))then
          target=log(max(m,tiny(1.0_dp)))-meanlog;a=max(shape(j),0.1_dp)
          do step=1,30
            newa=max(1.0e-5_dp,a-(log(a)-digamma_approx(a)-target)/(1.0_dp/a-trigamma_approx(a)))
            if(abs(newa-a)<1.0e-10_dp*(1.0_dp+a))exit
            a=newa
          end do
          shape(j)=newa
        end if
        scale(j)=m/shape(j)
      end do
      if(cs.and..not.present(fixed_shape))then;common=sum(shape*nk)/sum(nk);shape=common;do j=1,k
        scale(j)=sum(post(:,j)*x)/max(nk(j)*common,tiny(1.0_dp));end do
      end if
      call estep(newll);if(status/=0)exit
      history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=lambda;result%shape=shape;result%scale=scale;result%posterior=post;result%loglik=ll
    result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations.and.status==0;result%status=merge(0,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      ll=0.0_dp;status=0
      do i=1,n
        do j=1,k;lw(j)=log(max(lambda(j),tiny(1.0_dp)))+gamma_logpdf(x(i),shape(j),scale(j));end do
        call normalize_logweights(lw,post(i,:),ln);ll=ll+ln
      end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine gammamix_em

  subroutine multmix_em(y,k,result,control)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::k
    type(multinomial_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    type(em_control)::ctl
    integer::n,p,i,j,iter
    real(dp),allocatable::lambda(:),theta(:,:),post(:,:),history(:),lw(:),nk(:)
    real(dp)::ll,newll,ln,diff
    ctl=em_control();if(present(control))ctl=control;n=size(y,1);p=size(y,2)
    if(any(y<0.0_dp).or.n<k)then;result%status=MIXTOOLS_INVALID_ARGUMENT;return;end if
    allocate(lambda(k),theta(k,p),post(n,k),history(ctl%max_iterations+1),lw(k),nk(k))
    lambda=1.0_dp/real(k,dp)
    do j=1,k;theta(j,:)=y(1+mod((j-1)*max(1,n/k),n),:)+1.0_dp;theta(j,:)=theta(j,:)/sum(theta(j,:));end do
    call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1);lambda=nk/real(n,dp)
      do j=1,k;theta(j,:)=matmul(post(:,j),y)+1.0e-12_dp;theta(j,:)=theta(j,:)/sum(theta(j,:));end do
      call estep(newll);history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=lambda;result%theta=theta;result%posterior=post;result%loglik=ll
    result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations;result%status=merge(0,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      ll=0.0_dp
      do i=1,n
        do j=1,k;lw(j)=log(max(lambda(j),tiny(1.0_dp)))+multinomial_logpmf(y(i,:),theta(j,:));end do
        call normalize_logweights(lw,post(i,:),ln);ll=ll+ln
      end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine multmix_em

  subroutine repnormmix_em(x,k,result,control,common_mean,common_sigma)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::k
    type(mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    logical,intent(in),optional::common_mean,common_sigma
    type(em_control)::ctl
    logical::cm,cs
    integer::m,n,i,j,iter,status
    real(dp),allocatable::lambda(:),mu(:),sigma(:),post(:,:),history(:),lw(:),nk(:)
    real(dp)::ll,newll,ln,diff,sumsq,gm,gv
    ctl=em_control();if(present(control))ctl=control;cm=.false.;cs=.false.
    if(present(common_mean))cm=common_mean;if(present(common_sigma))cs=common_sigma
    m=size(x,1);n=size(x,2)
    call initial_univariate(reshape(x,[m*n]),k,lambda,mu,sigma,status)
    allocate(post(n,k),history(ctl%max_iterations+1),lw(k),nk(k));call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1);lambda=nk/real(n,dp)
      if(cm)then;gm=sum(x)/real(m*n,dp);mu=gm;else
        do j=1,k;mu(j)=sum(spread(post(:,j),1,m)*x)/max(real(m,dp)*nk(j),tiny(1.0_dp));end do
      end if
      if(cs)then
        gv=0.0_dp;do j=1,k;gv=gv+sum(spread(post(:,j),1,m)*(x-mu(j))**2);end do
        sigma=sqrt(max(gv/real(m*n,dp),ctl%minimum_scale**2))
      else
        do j=1,k;sigma(j)=sqrt(max(sum(spread(post(:,j),1,m)*(x-mu(j))**2) &
          /max(real(m,dp)*nk(j),tiny(1.0_dp)),ctl%minimum_scale**2));end do
      end if
      call estep(newll);history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=lambda;result%mu=mu;result%sigma=sigma;result%posterior=post;result%loglik=ll
    result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations;result%status=merge(0,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      ll=0.0_dp
      do i=1,n
        do j=1,k
          sumsq=sum((x(:,i)-mu(j))**2)
          lw(j)=log(max(lambda(j),tiny(1.0_dp)))-real(m,dp)*log(sigma(j))-0.5_dp*sumsq/sigma(j)**2
        end do
        call normalize_logweights(lw,post(i,:),ln);ll=ll+ln-0.5_dp*real(m,dp)*log(2.0_dp*acos(-1.0_dp))
      end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine repnormmix_em
end module mixtools_parametric
