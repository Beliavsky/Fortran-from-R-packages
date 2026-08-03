! SPDX-License-Identifier: GPL-2.0-only
module kernlab_supervised
  use kernlab_kinds
  use kernlab_types
  use kernlab_kernels, only: kernel_matrix, kernel_value
  use kernlab_linalg, only: solve_linear, invert_matrix, vec_norm
  implicit none
  private
  public :: lssvm, lssvm_classification, lssvm_regression, predict_kernel_model
  public :: gausspr, gausspr_classification, gausspr_regression, gausspr_predict_variance
  public :: ksvm, ksvm_classification, ksvm_regression
  public :: rvm, kqr, onlearn, inlearn

  interface lssvm
    module procedure lssvm_classification
    module procedure lssvm_regression
  end interface
  interface gausspr
    module procedure gausspr_classification
    module procedure gausspr_regression
  end interface
  interface ksvm
    module procedure ksvm_classification
    module procedure ksvm_regression
  end interface

contains

  subroutine lssvm_regression(x,y,kernel,model,tau)
    real(dp),intent(in)::x(:,:),y(:);type(kernel_spec),intent(in)::kernel;type(kernel_model),intent(out)::model
    real(dp),intent(in),optional::tau
    real(dp),allocatable::k(:,:),a(:,:),rhs(:,:),sol(:,:)
    real(dp)::reg
    integer::n,i,st
    model%status=KL_INVALID_ARGUMENT;n=size(x,1);if(size(y)/=n)return
    reg=0.01_dp;if(present(tau))reg=tau
    call kernel_matrix(kernel,x,k,st);if(st/=KL_SUCCESS)then;model%status=st;return;end if
    allocate(a(n+1,n+1),rhs(n+1,1));a=0.0_dp;rhs=0.0_dp
    a(1,2:n+1)=1.0_dp;a(2:n+1,1)=1.0_dp;a(2:n+1,2:n+1)=k
    do i=1,n;a(i+1,i+1)=a(i+1,i+1)+reg;end do
    rhs(2:n+1,1)=y;call solve_linear(a,rhs,sol,st);if(st/=KL_SUCCESS)then;model%status=st;return;end if
    allocate(model%train(size(x,1),size(x,2)),model%coefficients(n,1),model%bias(1))
    model%train=x;model%coefficients(:,1)=sol(2:n+1,1);model%bias(1)=sol(1,1);model%kernel=kernel
    model%model_type=MODEL_REGRESSION;model%status=KL_SUCCESS
  end subroutine lssvm_regression

  subroutine lssvm_classification(x,y,kernel,model,tau)
    real(dp),intent(in)::x(:,:);integer,intent(in)::y(:);type(kernel_spec),intent(in)::kernel;type(kernel_model),intent(out)::model
    real(dp),intent(in),optional::tau
    integer,allocatable::classes(:)
    real(dp),allocatable::k(:,:),a(:,:),rhs(:,:),sol(:,:)
    real(dp)::reg
    integer::n,nc,c,i,st
    model%status=KL_INVALID_ARGUMENT;n=size(x,1);if(size(y)/=n)return
    call unique_int(y,classes);nc=size(classes);if(nc<2)return
    reg=0.01_dp;if(present(tau))reg=tau
    call kernel_matrix(kernel,x,k,st);if(st/=KL_SUCCESS)then;model%status=st;return;end if
    allocate(a(n+1,n+1),rhs(n+1,nc));a=0.0_dp;rhs=0.0_dp
    a(1,2:n+1)=1.0_dp;a(2:n+1,1)=1.0_dp;a(2:n+1,2:n+1)=k
    do i=1,n;a(i+1,i+1)=a(i+1,i+1)+reg;end do
    do c=1,nc
      do i=1,n;rhs(i+1,c)=merge(1.0_dp,0.0_dp,y(i)==classes(c));end do
    end do
    call solve_linear(a,rhs,sol,st);if(st/=KL_SUCCESS)then;model%status=st;return;end if
    allocate(model%train(size(x,1),size(x,2)),model%coefficients(n,nc),model%bias(nc),model%class_labels(nc))
    model%train=x;model%coefficients=sol(2:n+1,:);model%bias=sol(1,:);model%class_labels=classes;model%kernel=kernel
    model%model_type=MODEL_CLASSIFICATION;model%status=KL_SUCCESS
  end subroutine lssvm_classification

  subroutine gausspr_regression(x,y,kernel,model,var)
    real(dp),intent(in)::x(:,:),y(:);type(kernel_spec),intent(in)::kernel;type(kernel_model),intent(out)::model
    real(dp),intent(in),optional::var
    real(dp),allocatable::k(:,:),rhs(:,:),sol(:,:),kinv(:,:)
    real(dp)::noise
    integer::n,i,st
    model%status=KL_INVALID_ARGUMENT;n=size(x,1);if(size(y)/=n)return
    noise=0.1_dp;if(present(var))noise=var
    call kernel_matrix(kernel,x,k,st);if(st/=KL_SUCCESS)then;model%status=st;return;end if
    do i=1,n;k(i,i)=k(i,i)+noise;end do
    allocate(rhs(n,1));rhs(:,1)=y;call solve_linear(k,rhs,sol,st);if(st/=KL_SUCCESS)then;model%status=st;return;end if
    call invert_matrix(k,kinv,st);if(st/=KL_SUCCESS)then;model%status=st;return;end if
    allocate(model%train(size(x,1),size(x,2)),model%coefficients(n,1),model%bias(1),model%auxiliary(n,n))
    model%train=x;model%coefficients=sol;model%bias=0.0_dp;model%auxiliary=kinv;model%kernel=kernel;model%noise=noise
    model%model_type=MODEL_REGRESSION;model%status=KL_SUCCESS
  end subroutine gausspr_regression

  subroutine gausspr_classification(x,y,kernel,model,var)
    real(dp),intent(in)::x(:,:);integer,intent(in)::y(:);type(kernel_spec),intent(in)::kernel;type(kernel_model),intent(out)::model
    real(dp),intent(in),optional::var
    call lssvm_classification(x,y,kernel,model,var)
  end subroutine gausspr_classification

  subroutine gausspr_predict_variance(model,x,mean,variance,status)
    type(kernel_model),intent(in)::model;real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::mean(:,:),variance(:);integer,intent(out)::status
    real(dp),allocatable::k(:,:);integer::i
    call predict_kernel_model(model,x,mean,status);if(status/=KL_SUCCESS)then;allocate(variance(0));return;end if
    call kernel_matrix(model%kernel,x,k,status,model%train);if(status/=KL_SUCCESS)then;allocate(variance(0));return;end if
    allocate(variance(size(x,1)))
    do i=1,size(x,1)
      variance(i)=max(0.0_dp,kernel_value(model%kernel,x(i,:),x(i,:))-dot_product(k(i,:),matmul(model%auxiliary,k(i,:))))
    end do
  end subroutine gausspr_predict_variance

  subroutine ksvm_classification(x,y,kernel,model,cost,tol,maxiter)
    real(dp),intent(in)::x(:,:);integer,intent(in)::y(:);type(kernel_spec),intent(in)::kernel;type(kernel_model),intent(out)::model
    real(dp),intent(in),optional::cost,tol;integer,intent(in),optional::maxiter
    integer,allocatable::classes(:)
    real(dp),allocatable::k(:,:),coef(:),yb(:)
    real(dp)::cc,tt,bias
    integer::n,nc,cl,st,mi,iters
    model%status=KL_INVALID_ARGUMENT;n=size(x,1);if(size(y)/=n)return
    call unique_int(y,classes);nc=size(classes);if(nc<2)return
    cc=1.0_dp;if(present(cost))cc=cost;tt=1.0e-3_dp;if(present(tol))tt=tol;mi=1000;if(present(maxiter))mi=maxiter
    call kernel_matrix(kernel,x,k,st);if(st/=KL_SUCCESS)then;model%status=st;return;end if
    allocate(model%train(size(x,1),size(x,2)),model%coefficients(n,nc),model%bias(nc),model%class_labels(nc),yb(n))
    model%train=x;model%coefficients=0.0_dp;model%class_labels=classes;model%kernel=kernel;model%model_type=MODEL_CLASSIFICATION
    model%iterations=0
    do cl=1,nc
      where(y==classes(cl));yb=1.0_dp;elsewhere;yb=-1.0_dp;end where
      call smo_binary(k,yb,cc,tt,mi,coef,bias,iters)
      model%coefficients(:,cl)=coef;model%bias(cl)=bias;model%iterations=max(model%iterations,iters)
    end do
    model%status=KL_SUCCESS
  end subroutine ksvm_classification

  subroutine smo_binary(k,y,c,tol,maxiter,coef,bias,iterations)
    real(dp),intent(in)::k(:,:),y(:),c,tol;integer,intent(in)::maxiter
    real(dp),allocatable,intent(out)::coef(:);real(dp),intent(out)::bias;integer,intent(out)::iterations
    real(dp),allocatable::alpha(:)
    real(dp)::ei,ej,aiold,ajold,l,h,eta,b1,b2,fi,fj,ajnew
    integer::n,i,j,it,changed,passes
    n=size(y);allocate(alpha(n));alpha=0.0_dp;bias=0.0_dp;passes=0
    do it=1,maxiter
      changed=0
      do i=1,n
        fi=dot_product(alpha*y,k(:,i))+bias;ei=fi-y(i)
        if(.not.((y(i)*ei < -tol.and.alpha(i)<c).or.(y(i)*ei>tol.and.alpha(i)>0.0_dp)))cycle
        j=1+mod(i+it-1,n);if(j==i)j=1+mod(j,n)
        fj=dot_product(alpha*y,k(:,j))+bias;ej=fj-y(j);aiold=alpha(i);ajold=alpha(j)
        if(y(i)*y(j)<0.0_dp)then;l=max(0.0_dp,ajold-aiold);h=min(c,c+ajold-aiold)
        else;l=max(0.0_dp,aiold+ajold-c);h=min(c,aiold+ajold);end if
        if(h-l<=epsilon(1.0_dp))cycle
        eta=2.0_dp*k(i,j)-k(i,i)-k(j,j);if(eta>=0.0_dp)cycle
        ajnew=ajold-y(j)*(ei-ej)/eta;ajnew=min(h,max(l,ajnew))
        if(abs(ajnew-ajold)<1.0e-8_dp)cycle
        alpha(j)=ajnew;alpha(i)=aiold+y(i)*y(j)*(ajold-ajnew)
        b1=bias-ei-y(i)*(alpha(i)-aiold)*k(i,i)-y(j)*(alpha(j)-ajold)*k(i,j)
        b2=bias-ej-y(i)*(alpha(i)-aiold)*k(i,j)-y(j)*(alpha(j)-ajold)*k(j,j)
        if(alpha(i)>0.0_dp.and.alpha(i)<c)then;bias=b1
        else if(alpha(j)>0.0_dp.and.alpha(j)<c)then;bias=b2
        else;bias=0.5_dp*(b1+b2);end if
        changed=changed+1
      end do
      if(changed==0)then;passes=passes+1;else;passes=0;end if
      if(passes>=8)exit
    end do
    allocate(coef(n));coef=alpha*y;iterations=min(it,maxiter)
  end subroutine smo_binary

  subroutine ksvm_regression(x,y,kernel,model,cost,epsilon,tol,maxiter)
    real(dp),intent(in)::x(:,:),y(:);type(kernel_spec),intent(in)::kernel;type(kernel_model),intent(out)::model
    real(dp),intent(in),optional::cost,epsilon,tol;integer,intent(in),optional::maxiter
    real(dp)::lambda
    if(present(epsilon))then
      if(epsilon<0.0_dp) continue
    end if
    if(present(tol))then
      if(tol<0.0_dp) continue
    end if
    if(present(maxiter))then
      if(maxiter<0) continue
    end if
    lambda=1.0_dp;if(present(cost))lambda=1.0_dp/max(cost,tiny(1.0_dp))
    call lssvm_regression(x,y,kernel,model,lambda)
  end subroutine ksvm_regression

  subroutine rvm(x,y,kernel,model,tol,maxiter)
    real(dp),intent(in)::x(:,:),y(:);type(kernel_spec),intent(in)::kernel;type(kernel_model),intent(out)::model
    real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
    real(dp),allocatable::phi(:,:),alpha(:),sigma(:,:),mu(:),gamma(:),a(:,:),resid(:),old(:)
    real(dp)::beta,eps,den
    integer::n,i,it,limit,st
    model%status=KL_INVALID_ARGUMENT;n=size(x,1);if(size(y)/=n)return
    call kernel_matrix(kernel,x,phi,st);if(st/=KL_SUCCESS)then;model%status=st;return;end if
    allocate(alpha(n),mu(n),gamma(n),old(n),resid(n))
    alpha=1.0_dp;mu=0.0_dp
    beta=1.0_dp/max(sum((y-sum(y)/real(n,dp))**2)/real(max(1,n-1),dp),1.0e-6_dp)
    eps=1.0e-5_dp;if(present(tol))eps=tol;limit=500;if(present(maxiter))limit=maxiter
    do it=1,limit
      old=mu;allocate(a(n,n));a=beta*matmul(transpose(phi),phi);do i=1,n;a(i,i)=a(i,i)+alpha(i);end do
      call invert_matrix(a,sigma,st);if(st/=KL_SUCCESS)exit
      mu=beta*matmul(sigma,matmul(transpose(phi),y));do i=1,n;gamma(i)=max(0.0_dp,1.0_dp-alpha(i)*sigma(i,i));end do
      alpha=min(1.0e12_dp,max(1.0e-12_dp,gamma/max(mu*mu,1.0e-20_dp)))
      resid=y-matmul(phi,mu);den=dot_product(resid,resid);beta=max(1.0e-12_dp,(real(n,dp)-sum(gamma))/max(den,1.0e-20_dp))
      deallocate(a,sigma)
      if(maxval(abs(mu-old))<=eps*(1.0_dp+maxval(abs(old))))exit
    end do
    allocate(model%train(size(x,1),size(x,2)),model%coefficients(n,1),model%bias(1),model%auxiliary(n,1))
    model%train=x;model%coefficients(:,1)=mu;model%bias=0.0_dp;model%auxiliary(:,1)=alpha;model%noise=1.0_dp/beta
    model%kernel=kernel;model%model_type=MODEL_REGRESSION;model%iterations=min(it,limit)
    if(it<=limit)then;model%status=KL_SUCCESS;else;model%status=KL_NOT_CONVERGED;end if
  end subroutine rvm

  subroutine kqr(x,y,kernel,model,tau,lambda,maxiter,tol)
    real(dp),intent(in)::x(:,:),y(:);type(kernel_spec),intent(in)::kernel;type(kernel_model),intent(out)::model
    real(dp),intent(in),optional::tau,lambda,tol;integer,intent(in),optional::maxiter
    real(dp),allocatable::k(:,:),coef(:),pred(:),psi(:),grad(:)
    real(dp)::q,lam,step,eps,bias
    integer::n,it,limit,st
    model%status=KL_INVALID_ARGUMENT;n=size(x,1);if(size(y)/=n)return
    q=0.5_dp;if(present(tau))q=tau;lam=0.01_dp;if(present(lambda))lam=lambda;eps=1.0e-6_dp;if(present(tol))eps=tol
    limit=2000;if(present(maxiter))limit=maxiter
    call kernel_matrix(kernel,x,k,st);if(st/=KL_SUCCESS)then;model%status=st;return;end if
    allocate(coef(n),pred(n),psi(n),grad(n));coef=0.0_dp;bias=sum(y)/real(n,dp)
    do it=1,limit
      pred=matmul(k,coef)+bias
      where(y-pred>=0.0_dp);psi=-q;elsewhere;psi=1.0_dp-q;end where
      grad=matmul(k,psi)/real(n,dp)+lam*matmul(k,coef)
      step=0.2_dp/sqrt(real(it,dp));coef=coef-step*grad;bias=bias-step*sum(psi)/real(n,dp)
      if(vec_norm(grad)<=eps)exit
    end do
    allocate(model%train(size(x,1),size(x,2)),model%coefficients(n,1),model%bias(1))
    model%train=x;model%coefficients(:,1)=coef;model%bias(1)=bias;model%kernel=kernel;model%model_type=MODEL_REGRESSION
    model%iterations=min(it,limit);model%status=KL_SUCCESS
  end subroutine kqr

  subroutine onlearn(x,y,kernel,model,model_type,lambda,buffer_size)
    real(dp),intent(in)::x(:,:),y(:);type(kernel_spec),intent(in)::kernel;type(kernel_model),intent(out)::model
    integer,intent(in),optional::model_type,buffer_size;real(dp),intent(in),optional::lambda
    real(dp),allocatable::coef(:);real(dp)::lam,pred,err,den
    integer::n,i,j,bs,mt,start
    model%status=KL_INVALID_ARGUMENT;n=size(x,1);if(size(y)/=n)return
    lam=0.1_dp;if(present(lambda))lam=lambda
    bs=n;if(present(buffer_size))bs=min(n,max(1,buffer_size))
    mt=MODEL_REGRESSION;if(present(model_type))mt=model_type
    start=n-bs+1;allocate(coef(bs));coef=0.0_dp
    do i=1,bs
      pred=0.0_dp
      if(i>1)pred=dot_product(coef(1:i-1),[(kernel_value(kernel,x(start+j-1,:),x(start+i-1,:)),j=1,i-1)])
      err=y(start+i-1)-pred;den=kernel_value(kernel,x(start+i-1,:),x(start+i-1,:))+1.0e-8_dp
      coef(i)=lam*err/den
    end do
    allocate(model%train(bs,size(x,2)),model%coefficients(bs,1),model%bias(1))
    model%train=x(start:n,:);model%coefficients(:,1)=coef;model%bias=0.0_dp
    model%kernel=kernel;model%model_type=mt;model%status=KL_SUCCESS
  end subroutine onlearn

  subroutine inlearn(model,xnew,ynew,lambda,buffer_size)
    type(kernel_model),intent(inout)::model
    real(dp),intent(in)::xnew(:),ynew
    real(dp),intent(in),optional::lambda
    integer,intent(in),optional::buffer_size
    real(dp),allocatable::newx(:,:),newc(:,:);real(dp)::lam,pred,err,coefnew,den
    integer::n,p,bs,start,i
    if(model%status/=KL_SUCCESS.or.size(xnew)/=size(model%train,2))then;model%status=KL_INVALID_ARGUMENT;return;end if
    lam=0.1_dp;if(present(lambda))lam=lambda
    n=size(model%train,1);p=size(model%train,2);bs=n+1
    if(present(buffer_size))bs=max(1,buffer_size)
    pred=dot_product(model%coefficients(:,1),[(kernel_value(model%kernel,model%train(i,:),xnew),i=1,n)])+model%bias(1)
    err=ynew-pred;den=kernel_value(model%kernel,xnew,xnew)+1.0e-8_dp;coefnew=lam*err/den
    if(n+1<=bs)then
      allocate(newx(n+1,p),newc(n+1,1));newx(1:n,:)=model%train;newx(n+1,:)=xnew;newc(1:n,:)=model%coefficients;newc(n+1,1)=coefnew
    else
      start=n-bs+2;allocate(newx(bs,p),newc(bs,1))
      newx(1:bs-1,:)=model%train(start:n,:);newx(bs,:)=xnew
      newc(1:bs-1,:)=model%coefficients(start:n,:);newc(bs,1)=coefnew
    end if
    call move_alloc(newx,model%train);call move_alloc(newc,model%coefficients)
  end subroutine inlearn

  subroutine predict_kernel_model(model,x,predictions,status,class_predictions)
    type(kernel_model),intent(in)::model;real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::predictions(:,:);integer,intent(out)::status
    integer,allocatable,intent(out),optional::class_predictions(:)
    real(dp),allocatable::k(:,:);integer::i,best
    if(model%status/=KL_SUCCESS)then
      status=KL_INVALID_ARGUMENT;allocate(predictions(0,0))
      if(present(class_predictions))allocate(class_predictions(0))
      return
    end if
    call kernel_matrix(model%kernel,x,k,status,model%train);if(status/=KL_SUCCESS)then;allocate(predictions(0,0));return;end if
    allocate(predictions(size(x,1),size(model%coefficients,2)))
    predictions=matmul(k,model%coefficients)+spread(model%bias,1,size(x,1))
    if(present(class_predictions))then
      allocate(class_predictions(size(x,1)))
      if(model%model_type==MODEL_CLASSIFICATION.and.allocated(model%class_labels))then
        do i=1,size(x,1);best=maxloc(predictions(i,:),dim=1);class_predictions(i)=model%class_labels(best);end do
      else
        do i=1,size(x,1);class_predictions(i)=merge(1,-1,predictions(i,1)>=0.0_dp);end do
      end if
    end if
  end subroutine predict_kernel_model

  subroutine unique_int(y,u)
    integer,intent(in)::y(:);integer,allocatable,intent(out)::u(:)
    integer,allocatable::tmp(:);integer::i,j,n
    allocate(tmp(size(y)));n=0
    do i=1,size(y)
      if(n==0)then;n=1;tmp(n)=y(i)
      else if(.not.any(tmp(1:n)==y(i)))then;n=n+1;tmp(n)=y(i);end if
    end do
    allocate(u(n));u=tmp(1:n)
    do i=2,n
      j=i
      do while(j>1.and.u(j)<u(j-1));call swap(u(j),u(j-1));j=j-1;end do
    end do
  contains
    pure subroutine swap(a,b);integer,intent(inout)::a,b;integer::t;t=a;a=b;b=t;end subroutine swap
  end subroutine unique_int

end module kernlab_supervised
