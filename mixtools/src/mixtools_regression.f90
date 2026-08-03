! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_regression
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mixtools_kinds, only : dp, pi
  use mixtools_status
  use mixtools_types
  use mixtools_distributions, only : normal_logpdf, normalize_logweights, logistic
  use mixtools_distributions, only : binomial_logpmf, poisson_logpmf
  use mixtools_linalg, only : solve_least_squares
  use mixtools_utilities, only : kernel_value
  implicit none
  private
  public :: regmix_em, regmix_em_lambda, regmix_em_loc, logisregmix_em
  public :: poisregmix_em, segregmix_em, hme_em, flaremix_em
  public :: regmix_em_mixed
contains
  subroutine add_intercept_matrix(x,addintercept,xa)
    real(dp),intent(in)::x(:,:)
    logical,intent(in)::addintercept
    real(dp),allocatable,intent(out)::xa(:,:)
    if(addintercept)then
      allocate(xa(size(x,1),size(x,2)+1));xa(:,1)=1.0_dp;xa(:,2:)=x
    else
      xa=x
    end if
  end subroutine add_intercept_matrix

  subroutine weighted_ls(x,y,w,beta,status,ridge)
    real(dp),intent(in)::x(:,:),y(:),w(:)
    real(dp),allocatable,intent(out)::beta(:)
    integer,intent(out)::status
    real(dp),intent(in),optional::ridge
    real(dp),allocatable::xw(:,:),yw(:),cov(:,:)
    real(dp)::rss,r
    integer::i
    if(size(x,1)/=size(y).or.size(w)/=size(y))then;allocate(beta(0));status=MIXTOOLS_DIMENSION_ERROR;return;end if
    allocate(xw(size(x,1),size(x,2)),yw(size(y)));do i=1,size(y)
      xw(i,:)=sqrt(max(w(i),0.0_dp))*x(i,:);yw(i)=sqrt(max(w(i),0.0_dp))*y(i)
    end do
    r=1.0e-8_dp;if(present(ridge))r=ridge
    if(r>0.0_dp)then
      call ridge_ls(xw,yw,r,beta,status)
    else
      call solve_least_squares(xw,yw,beta,cov,rss,status)
    end if
  contains
    subroutine ridge_ls(a,b,lambda,coef,stat)
      real(dp),intent(in)::a(:,:),b(:),lambda
      real(dp),allocatable,intent(out)::coef(:)
      integer,intent(out)::stat
      real(dp),allocatable::aa(:,:),bb(:),cv(:,:)
      real(dp)::rr
      integer::n,p,ii
      n=size(a,1);p=size(a,2);allocate(aa(n+p,p),bb(n+p));aa(1:n,:)=a;bb(1:n)=b
      aa(n+1:,:)=0.0_dp;bb(n+1:)=0.0_dp
      do ii=1,p;aa(n+ii,ii)=sqrt(lambda);end do
      call solve_least_squares(aa,bb,coef,cv,rr,stat)
    end subroutine ridge_ls
  end subroutine weighted_ls

  subroutine initialize_regression(y,x,k,beta,sigma,lambda,status)
    real(dp),intent(in)::y(:),x(:,:)
    integer,intent(in)::k
    real(dp),allocatable,intent(out)::beta(:,:),sigma(:),lambda(:)
    integer,intent(out)::status
    real(dp),allocatable::b(:),w(:)
    real(dp)::sd
    integer::j
    allocate(beta(size(x,2),k),sigma(k),lambda(k),w(size(y)));w=1.0_dp
    call weighted_ls(x,y,w,b,status)
    if(status/=MIXTOOLS_SUCCESS)return
    sd=sqrt(sum((y-matmul(x,b))**2)/real(max(1,size(y)-size(x,2)),dp));sd=max(sd,1.0e-3_dp)
    do j=1,k
      beta(:,j)=b;beta(1,j)=beta(1,j)+(real(j,dp)-(real(k,dp)+1.0_dp)/2.0_dp)*0.5_dp*sd
      sigma(j)=sd;lambda(j)=1.0_dp/real(k,dp)
    end do
  end subroutine initialize_regression

  subroutine regmix_em(y,x,k,result,control,addintercept,common_beta,common_sigma)
    real(dp),intent(in)::y(:),x(:,:)
    integer,intent(in)::k
    type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    logical,intent(in),optional::addintercept,common_beta,common_sigma
    type(em_control)::ctl
    logical::ai,cb,cs
    real(dp),allocatable::xa(:,:),beta(:,:),sigma(:),lambda(:),post(:,:),history(:),lw(:),nk(:),b(:),w(:)
    real(dp)::ll,newll,ln,diff,mu,pooled
    integer::n,p,i,j,iter,status
    ctl=em_control();if(present(control))ctl=control;ai=.true.;cb=.false.;cs=.false.
    if(present(addintercept))ai=addintercept;if(present(common_beta))cb=common_beta;if(present(common_sigma))cs=common_sigma
    call add_intercept_matrix(x,ai,xa);n=size(y);p=size(xa,2)
    if(size(xa,1)/=n.or.n<k)then;result%status=MIXTOOLS_DIMENSION_ERROR;return;end if
    call initialize_regression(y,xa,k,beta,sigma,lambda,status);if(status/=0)then;result%status=status;return;end if
    allocate(post(n,k),history(ctl%max_iterations+1),lw(k),nk(k),w(n));call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1);lambda=max(nk/real(n,dp),tiny(1.0_dp));lambda=lambda/sum(lambda)
      if(cb)then
        w=sum(post,dim=2);call weighted_ls(xa,y,w,b,status,ctl%ridge);if(status/=0)exit
        beta=spread(b,2,k)
      else
        do j=1,k;call weighted_ls(xa,y,post(:,j),b,status,ctl%ridge);if(status/=0)exit;beta(:,j)=b;end do
        if(status/=0)exit
      end if
      if(cs)then
        pooled=0.0_dp;do j=1,k;pooled=pooled+sum(post(:,j)*(y-matmul(xa,beta(:,j)))**2);end do
        sigma=sqrt(max(pooled/real(n,dp),ctl%minimum_scale**2))
      else
        do j=1,k;sigma(j)=sqrt(max(sum(post(:,j)*(y-matmul(xa,beta(:,j)))**2)/max(nk(j),tiny(1.0_dp)), &
          ctl%minimum_scale**2));end do
      end if
      call estep(newll);if(status/=0)exit;history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=lambda;result%beta=beta;result%sigma=sigma;result%posterior=post;result%loglik=ll
    result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations.and.status==0;result%status=merge(0,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      ll=0.0_dp;status=0
      do i=1,n
        do j=1,k
          mu=dot_product(xa(i,:),beta(:,j));lw(j)=log(max(lambda(j),tiny(1.0_dp)))+normal_logpdf(y(i),mu,sigma(j))
        end do
        call normalize_logweights(lw,post(i,:),ln);if(.not.ieee_is_finite(ln))then;status=MIXTOOLS_NUMERICAL_ERROR;return;end if
        ll=ll+ln
      end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine regmix_em

  subroutine regmix_em_lambda(y,x,lambda_x,result,control,addintercept)
    real(dp),intent(in)::y(:),x(:,:),lambda_x(:,:)
    type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    logical,intent(in),optional::addintercept
    type(em_control)::ctl
    logical::ai
    real(dp),allocatable::xa(:,:),beta(:,:),sigma(:),post(:,:),history(:),lw(:),nk(:),b(:)
    real(dp)::ll,newll,ln,diff,mu
    integer::n,k,i,j,iter,status
    ctl=em_control();if(present(control))ctl=control;ai=.true.;if(present(addintercept))ai=addintercept
    call add_intercept_matrix(x,ai,xa);n=size(y);k=size(lambda_x,2)
    if(size(lambda_x,1)/=n.or.any(lambda_x<0.0_dp))then;result%status=MIXTOOLS_DIMENSION_ERROR;return;end if
    call initialize_regression(y,xa,k,beta,sigma,result%lambda,status);if(status/=0)then;result%status=status;return;end if
    deallocate(result%lambda);allocate(post(n,k),history(ctl%max_iterations+1),lw(k),nk(k));call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1)
      do j=1,k;call weighted_ls(xa,y,post(:,j),b,status,ctl%ridge);if(status/=0)exit;beta(:,j)=b
        sigma(j)=sqrt(max(sum(post(:,j)*(y-matmul(xa,beta(:,j)))**2)/max(nk(j),tiny(1.0_dp)),ctl%minimum_scale**2))
      end do
      if(status/=0)exit;call estep(newll);history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=sum(lambda_x,dim=1)/real(n,dp);result%beta=beta;result%sigma=sigma;result%posterior=post
    result%auxiliary=lambda_x;result%loglik=ll;result%iterations=min(iter,ctl%max_iterations)
    result%loglik_history=history(:result%iterations+1);result%converged=iter<=ctl%max_iterations.and.status==0
    result%status=merge(0,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      ll=0.0_dp;status=0
      do i=1,n
        do j=1,k;mu=dot_product(xa(i,:),beta(:,j));lw(j)=log(max(lambda_x(i,j),tiny(1.0_dp))) &
          +normal_logpdf(y(i),mu,sigma(j));end do
        call normalize_logweights(lw,post(i,:),ln);ll=ll+ln
      end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine regmix_em_lambda

  subroutine regmix_em_loc(y,x,k,bandwidth,result,control,addintercept,kernel)
    real(dp),intent(in)::y(:),x(:,:)
    integer,intent(in)::k
    real(dp),intent(in)::bandwidth
    type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    logical,intent(in),optional::addintercept
    integer,intent(in),optional::kernel
    type(regression_mixture_result)::base,current
    real(dp),allocatable::lambda_x(:,:)
    integer::i,j,h,n,kk
    real(dp)::den,w
    kk=1;if(present(kernel))kk=kernel;n=size(y)
    call regmix_em(y,x,k,base,control,addintercept)
    if(base%status/=0.and.base%status/=MIXTOOLS_NOT_CONVERGED)then;result=base;return;end if
    allocate(lambda_x(n,k))
    do i=1,n;do j=1,k;den=0.0_dp;lambda_x(i,j)=0.0_dp
      do h=1,n;w=kernel_value(sqrt(sum((x(i,:)-x(h,:))**2))/bandwidth,kk);den=den+w
        lambda_x(i,j)=lambda_x(i,j)+w*base%posterior(h,j)
      end do
      lambda_x(i,j)=lambda_x(i,j)/max(den,tiny(1.0_dp))
    end do;lambda_x(i,:)=lambda_x(i,:)/sum(lambda_x(i,:));end do
    call regmix_em_lambda(y,x,lambda_x,current,control,addintercept);result=current
  end subroutine regmix_em_loc

  subroutine glm_component_irls(y,ntrials,x,w,beta,family,control,status)
    real(dp),intent(in)::y(:),ntrials(:),x(:,:),w(:)
    real(dp),intent(inout)::beta(:)
    integer,intent(in)::family
    type(em_control),intent(in)::control
    integer,intent(out)::status
    real(dp),allocatable::eta(:),mu(:),varw(:),z(:),bnew(:)
    integer::iter
    allocate(eta(size(y)),mu(size(y)),varw(size(y)),z(size(y)))
    do iter=1,100
      eta=matmul(x,beta)
      if(family==1)then
        mu=ntrials*logistic(eta);varw=w*max(mu*(1.0_dp-mu/max(ntrials,1.0_dp)),1.0e-8_dp)
        z=eta+(y-mu)/max(mu*(1.0_dp-mu/max(ntrials,1.0_dp)),1.0e-8_dp)
      else
        mu=exp(min(eta,30.0_dp));varw=w*max(mu,1.0e-8_dp);z=eta+(y-mu)/max(mu,1.0e-8_dp)
      end if
      call weighted_ls(x,z,varw,bnew,status,control%ridge);if(status/=0)return
      if(maxval(abs(bnew-beta))<control%tolerance*(1.0_dp+maxval(abs(beta))))then;beta=bnew;return;end if
      beta=bnew
    end do
    status=MIXTOOLS_SUCCESS
  end subroutine glm_component_irls

  subroutine glm_mix_em(y,x,k,result,family,ntrials,control,addintercept)
    real(dp),intent(in)::y(:),x(:,:)
    integer,intent(in)::k,family
    type(regression_mixture_result),intent(out)::result
    real(dp),intent(in),optional::ntrials(:)
    type(em_control),intent(in),optional::control
    logical,intent(in),optional::addintercept
    type(em_control)::ctl
    logical::ai
    real(dp),allocatable::xa(:,:),trials(:),lambda(:),beta(:,:),post(:,:),history(:),lw(:),nk(:)
    real(dp)::ll,newll,ln,diff,eta,pr
    integer::n,p,i,j,iter,status
    ctl=em_control();if(present(control))ctl=control;ai=.true.;if(present(addintercept))ai=addintercept
    call add_intercept_matrix(x,ai,xa);n=size(y);p=size(xa,2);allocate(trials(n));trials=1.0_dp
    if(present(ntrials))then;if(size(ntrials)/=n)then;result%status=MIXTOOLS_DIMENSION_ERROR;return;end if;trials=ntrials;end if
    if(family==1.and.any(y>trials))then;result%status=MIXTOOLS_INVALID_ARGUMENT;return;end if
    allocate(lambda(k),beta(p,k),post(n,k),history(ctl%max_iterations+1),lw(k),nk(k))
    lambda=1.0_dp/real(k,dp);beta=0.0_dp
    if(family==1)then
      do j=1,k
        beta(1,j)=log((sum(y)+real(j,dp))/(sum(trials-y)+real(k-j+1,dp)))
      end do
    else
      do j=1,k
        beta(1,j)=log(max((sum(y)+real(j,dp))/real(n+k,dp),1.0e-8_dp))
      end do
    end if
    call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1);lambda=max(nk/real(n,dp),tiny(1.0_dp));lambda=lambda/sum(lambda)
      do j=1,k;call glm_component_irls(y,trials,xa,post(:,j),beta(:,j),family,ctl,status);if(status/=0)exit;end do
      if(status/=0)exit;call estep(newll);history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=lambda;result%beta=beta;allocate(result%sigma(k));result%sigma=1.0_dp;result%posterior=post
    result%loglik=ll;result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations.and.status==0;result%status=merge(0,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      ll=0.0_dp;status=0
      do i=1,n;do j=1,k;eta=dot_product(xa(i,:),beta(:,j))
        if(family==1)then;pr=logistic(eta);lw(j)=log(max(lambda(j),tiny(1.0_dp)))+binomial_logpmf(y(i),trials(i),pr)
        else;lw(j)=log(max(lambda(j),tiny(1.0_dp)))+poisson_logpmf(y(i),exp(min(eta,30.0_dp)));end if
      end do;call normalize_logweights(lw,post(i,:),ln);ll=ll+ln;end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine glm_mix_em

  subroutine logisregmix_em(y,x,k,result,ntrials,control,addintercept)
    real(dp),intent(in)::y(:),x(:,:);integer,intent(in)::k
    type(regression_mixture_result),intent(out)::result
    real(dp),intent(in),optional::ntrials(:);type(em_control),intent(in),optional::control
    logical,intent(in),optional::addintercept
    call glm_mix_em(y,x,k,result,1,ntrials,control,addintercept)
  end subroutine logisregmix_em

  subroutine poisregmix_em(y,x,k,result,control,addintercept)
    real(dp),intent(in)::y(:),x(:,:);integer,intent(in)::k
    type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control;logical,intent(in),optional::addintercept
    call glm_mix_em(y,x,k,result,2,control=control,addintercept=addintercept)
  end subroutine poisregmix_em

  subroutine segregmix_em(y,x,segmented_column,k,breakpoints,result,control,addintercept)
    real(dp),intent(in)::y(:),x(:,:)
    integer,intent(in)::segmented_column,k
    real(dp),intent(in)::breakpoints(:)
    type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    logical,intent(in),optional::addintercept
    real(dp),allocatable::xaug(:,:),psi(:)
    type(regression_mixture_result)::fit,best
    real(dp)::lo,hi,cand,bestll
    integer::j,g,n,p
    if(size(breakpoints)/=k.or.segmented_column<1.or.segmented_column>size(x,2))then
      result%status=MIXTOOLS_DIMENSION_ERROR;return
    end if
    n=size(x,1);p=size(x,2);allocate(psi(k));psi=breakpoints;bestll=-huge(1.0_dp)
    do g=0,20
      allocate(xaug(n,p+k));xaug(:,1:p)=x
      do j=1,k;xaug(:,p+j)=max(0.0_dp,x(:,segmented_column)-psi(j));end do
      call regmix_em(y,xaug,k,fit,control,addintercept)
      if(fit%loglik>bestll)then;best=fit;bestll=fit%loglik;end if
      deallocate(xaug)
      if(g<20)then
        lo=minval(x(:,segmented_column));hi=maxval(x(:,segmented_column))
        do j=1,k;cand=lo+(hi-lo)*(real(g+1,dp)/21.0_dp);psi(j)=0.8_dp*psi(j)+0.2_dp*cand;end do
      end if
    end do
    result=best;allocate(result%auxiliary(1,k));result%auxiliary(1,:)=psi
  end subroutine segregmix_em

  subroutine hme_em(y,x,result,control,addintercept)
    real(dp),intent(in)::y(:),x(:,:)
    type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    logical,intent(in),optional::addintercept
    type(em_control)::ctl
    type(regression_mixture_result)::base,cur
    real(dp),allocatable::xa(:,:),wgate(:),lambda_x(:,:),grad(:),hess(:,:),delta(:),cov(:,:)
    real(dp)::rss,p1
    integer::i,iter,status,n,p
    ctl=em_control();if(present(control))ctl=control
    call regmix_em(y,x,2,base,ctl,addintercept);call add_intercept_matrix(x,.true.,xa);n=size(y);p=size(xa,2)
    allocate(wgate(p),lambda_x(n,2),grad(p),hess(p,p));wgate=0.0_dp
    do iter=1,min(100,ctl%max_iterations)
      do i=1,n;p1=logistic(dot_product(xa(i,:),wgate));lambda_x(i,:)=[p1,1.0_dp-p1];end do
      call regmix_em_lambda(y,x,lambda_x,cur,ctl,addintercept)
      grad=matmul(transpose(xa),cur%posterior(:,1)-lambda_x(:,1));hess=0.0_dp
      do i=1,n;hess=hess+lambda_x(i,1)*(1.0_dp-lambda_x(i,1))*outer(xa(i,:));end do
      call solve_least_squares(hess,grad,delta,cov,rss,status);if(status/=0)exit
      wgate=wgate+delta;if(maxval(abs(delta))<ctl%tolerance)exit
    end do
    result=cur
    if(allocated(result%auxiliary))deallocate(result%auxiliary)
    allocate(result%auxiliary(p,1));result%auxiliary(:,1)=wgate
  contains
    pure function outer(v) result(a)
      real(dp),intent(in)::v(:);real(dp)::a(size(v),size(v));integer::ii,jj
      do ii=1,size(v);do jj=1,size(v);a(ii,jj)=v(ii)*v(jj);end do;end do
    end function outer
  end subroutine hme_em

  subroutine flaremix_em(y,x,k,result,control,nu)
    real(dp),intent(in)::y(:),x(:,:)
    integer,intent(in)::k
    type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::nu
    type(em_control)::ctl
    real(dp)::barrier
    ctl=em_control();if(present(control))ctl=control;barrier=1.0_dp;if(present(nu))barrier=nu
    call regmix_em(y,x,k,result,ctl,.true.)
    if(allocated(result%sigma))result%sigma=max(result%sigma,barrier*ctl%minimum_scale)
  end subroutine flaremix_em

  subroutine regmix_em_mixed(y,x,groups,k,result,random_effects,control,addintercept)
    real(dp),intent(in)::y(:),x(:,:)
    integer,intent(in)::groups(:),k
    type(regression_mixture_result),intent(out)::result
    real(dp),allocatable,intent(out)::random_effects(:,:)
    type(em_control),intent(in),optional::control
    logical,intent(in),optional::addintercept
    type(em_control)::ctl
    logical::ai
    real(dp),allocatable::xf(:,:),design(:,:),coef(:,:),sigma(:),lambda(:),post(:,:),history(:),lw(:),nk(:),b(:)
    real(dp)::ll,newll,ln,diff,mu
    integer::n,p,ng,i,j,iter,status
    ctl=em_control();if(present(control))ctl=control;ai=.true.;if(present(addintercept))ai=addintercept
    if(size(groups)/=size(y).or.minval(groups)<1)then
      result%status=MIXTOOLS_DIMENSION_ERROR;allocate(random_effects(0,0));return
    end if
    call add_intercept_matrix(x,ai,xf);n=size(y);p=size(xf,2);ng=maxval(groups)
    allocate(design(n,p+ng));design(:,1:p)=xf;design(:,p+1:)=0.0_dp
    do i=1,n;design(i,p+groups(i))=1.0_dp;end do
    call initialize_regression(y,design,k,coef,sigma,lambda,status)
    if(status/=0)then;result%status=status;allocate(random_effects(0,0));return;end if
    allocate(post(n,k),history(ctl%max_iterations+1),lw(k),nk(k));call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1);lambda=max(nk/real(n,dp),tiny(1.0_dp));lambda=lambda/sum(lambda)
      do j=1,k
        call weighted_ls(design,y,post(:,j),b,status,max(ctl%ridge,1.0e-4_dp))
        if(status/=0)exit
        coef(:,j)=b
        sigma(j)=sqrt(max(sum(post(:,j)*(y-matmul(design,coef(:,j)))**2)/max(nk(j),tiny(1.0_dp)), &
          ctl%minimum_scale**2))
      end do
      if(status/=0)exit
      ! Center group effects within each component to separate them from the intercept.
      do j=1,k
        mu=sum(coef(p+1:p+ng,j))/real(ng,dp)
        coef(p+1:p+ng,j)=coef(p+1:p+ng,j)-mu
        if(ai)coef(1,j)=coef(1,j)+mu
      end do
      call estep(newll);history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=lambda;result%beta=coef(1:p,:);result%sigma=sigma;result%posterior=post
    allocate(random_effects(ng,k));random_effects=coef(p+1:p+ng,:);result%auxiliary=random_effects
    result%loglik=ll;result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations.and.status==0
    result%status=merge(MIXTOOLS_SUCCESS,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      ll=0.0_dp;status=0
      do i=1,n
        do j=1,k
          mu=dot_product(design(i,:),coef(:,j))
          lw(j)=log(max(lambda(j),tiny(1.0_dp)))+normal_logpdf(y(i),mu,sigma(j))
        end do
        call normalize_logweights(lw,post(i,:),ln);ll=ll+ln
      end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine regmix_em_mixed

end module mixtools_regression
