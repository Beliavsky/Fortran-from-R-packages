module mlr_learners
  use mlr_kinds, only : dp
  use mlr_types, only : linear_model, logistic_model, kmeans_model, knn_model
  use mlr_utils, only : mean_dp, median_dp, solve_linear, sigmoid
  use mlr_rng, only : rng_state, rng_integer
  implicit none
  private
  public :: fit_featureless_regression, predict_featureless_regression
  public :: fit_featureless_classifier, predict_featureless_classifier
  public :: fit_linear_regression, predict_linear_regression
  public :: fit_logistic_regression, predict_logistic_probability, predict_logistic_class
  public :: fit_kmeans, predict_kmeans
  public :: fit_knn_regression, fit_knn_classifier, predict_knn_regression, predict_knn_classifier
contains
  subroutine fit_featureless_regression(y, response, method)
    real(dp),intent(in)::y(:); real(dp),intent(out)::response; character(len=*),intent(in),optional::method
    character(len=16)::m
    m='mean';if(present(method))m=trim(method)
    if(trim(m)=='median')then;response=median_dp(y);else;response=mean_dp(y);end if
  end subroutine
  subroutine predict_featureless_regression(response,n,pred)
    real(dp),intent(in)::response;integer,intent(in)::n;real(dp),allocatable,intent(out)::pred(:)
    allocate(pred(n));pred=response
  end subroutine

  subroutine fit_featureless_classifier(y,nclass,probs)
    integer,intent(in)::y(:),nclass;real(dp),allocatable,intent(out)::probs(:)
    integer::k;allocate(probs(nclass));
    do k=1,nclass;probs(k)=real(count(y==k),dp)/real(size(y),dp);end do
  end subroutine
  subroutine predict_featureless_classifier(probs,n,response,prob)
    real(dp),intent(in)::probs(:);integer,intent(in)::n
    integer,allocatable,intent(out)::response(:);real(dp),allocatable,intent(out),optional::prob(:,:)
    integer::best,i
    best=maxloc(probs,dim=1);allocate(response(n));response=best
    if(present(prob))then;allocate(prob(n,size(probs)));do i=1,n;prob(i,:)=probs;end do;end if
  end subroutine

  subroutine fit_linear_regression(x,y,model,weights,ridge)
    real(dp),intent(in)::x(:,:),y(:);type(linear_model),intent(out)::model
    real(dp),intent(in),optional::weights(:),ridge
    real(dp),allocatable::z(:,:),a(:,:),b(:),sol(:),res(:);real(dp)::w,lam
    integer::i,j,k,p,info,n
    n=size(x,1);p=size(x,2)+1
    if(size(y)/=n)error stop 'fit_linear_regression: size mismatch'
    lam=0.0_dp;if(present(ridge))lam=max(0.0_dp,ridge)
    allocate(z(n,p));z(:,1)=1.0_dp;z(:,2:p)=x
    allocate(a(p,p),b(p));a=0.0_dp;b=0.0_dp
    do i=1,n
      w=1.0_dp;if(present(weights))w=weights(i)
      do j=1,p
        b(j)=b(j)+w*z(i,j)*y(i)
        do k=1,p;a(j,k)=a(j,k)+w*z(i,j)*z(i,k);end do
      end do
    end do
    do j=2,p;a(j,j)=a(j,j)+lam;end do
    call solve_linear(a,b,sol,info)
    if(info/=0)then
      do j=1,p;a(j,j)=a(j,j)+1.0e-10_dp;end do
      call solve_linear(a,b,sol,info)
    end if
    if(info/=0)error stop 'fit_linear_regression: singular system'
    model%coef=sol;allocate(res(n));res=y-matmul(z,sol)
    if(n>p)then;model%sigma=sqrt(sum(res*res)/real(n-p,dp));else;model%sigma=0.0_dp;end if
    model%fitted=.true.
  end subroutine
  subroutine predict_linear_regression(model,x,pred,se)
    type(linear_model),intent(in)::model;real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::pred(:);real(dp),allocatable,intent(out),optional::se(:)
    if(.not.model%fitted)error stop 'predict_linear_regression: model not fitted'
    allocate(pred(size(x,1)));pred=model%coef(1)+matmul(x,model%coef(2:))
    if(present(se))then;allocate(se(size(x,1)));se=model%sigma;end if
  end subroutine

  subroutine fit_logistic_regression(x,y,model,weights,maxiter,tol,ridge)
    real(dp),intent(in)::x(:,:);integer,intent(in)::y(:);type(logistic_model),intent(out)::model
    real(dp),intent(in),optional::weights(:),tol,ridge;integer,intent(in),optional::maxiter
    real(dp),allocatable::z(:,:),beta(:),a(:,:),b(:),step(:),pvec(:)
    real(dp)::w,lam,eps;integer::n,q,i,j,k,it,mi,info
    n=size(x,1);q=size(x,2)+1;if(size(y)/=n)error stop 'fit_logistic_regression: size mismatch'
    mi=100;if(present(maxiter))mi=maxiter;eps=1.0e-8_dp;if(present(tol))eps=tol
    lam=1.0e-8_dp;if(present(ridge))lam=max(0.0_dp,ridge)
    allocate(z(n,q));z(:,1)=1.0_dp;z(:,2:q)=x;allocate(beta(q));beta=0.0_dp
    allocate(a(q,q),b(q),pvec(n))
    model%converged=.false.
    do it=1,mi
      pvec=sigmoid(matmul(z,beta));a=0.0_dp;b=0.0_dp
      do i=1,n
        w=1.0_dp;if(present(weights))w=weights(i)
        do j=1,q
          b(j)=b(j)+w*z(i,j)*(real(y(i),dp)-pvec(i))
          do k=1,q;a(j,k)=a(j,k)+w*pvec(i)*(1.0_dp-pvec(i))*z(i,j)*z(i,k);end do
        end do
      end do
      do j=2,q;a(j,j)=a(j,j)+lam;end do
      call solve_linear(a,b,step,info);if(info/=0)exit
      beta=beta+step
      if(maxval(abs(step))<eps)then;model%converged=.true.;exit;end if
    end do
    model%coef=beta;model%iterations=min(it,mi)
  end subroutine
  subroutine predict_logistic_probability(model,x,prob)
    type(logistic_model),intent(in)::model;real(dp),intent(in)::x(:,:);real(dp),allocatable,intent(out)::prob(:)
    allocate(prob(size(x,1)));prob=sigmoid(model%coef(1)+matmul(x,model%coef(2:)))
  end subroutine
  subroutine predict_logistic_class(model,x,response,threshold)
    type(logistic_model),intent(in)::model;real(dp),intent(in)::x(:,:);integer,allocatable,intent(out)::response(:)
    real(dp),intent(in),optional::threshold;real(dp),allocatable::p(:);real(dp)::t;integer::i
    t=0.5_dp;if(present(threshold))t=threshold;call predict_logistic_probability(model,x,p)
    allocate(response(size(p)));do i=1,size(p);response(i)=merge(1,0,p(i)>=t);end do
  end subroutine

  subroutine fit_kmeans(x,k,model,rng,maxiter,tol)
    real(dp),intent(in)::x(:,:);integer,intent(in)::k;type(kmeans_model),intent(out)::model
    type(rng_state),intent(inout)::rng;integer,intent(in),optional::maxiter;real(dp),intent(in),optional::tol
    integer::n,p,i,c,it,mi,best;real(dp)::eps,d,bestd,shift
    integer,allocatable::assign(:),countc(:),chosen(:);real(dp),allocatable::newc(:,:)
    n=size(x,1);p=size(x,2);if(k<1.or.k>n)error stop 'fit_kmeans: invalid k'
    mi=100;if(present(maxiter))mi=maxiter;eps=1.0e-8_dp;if(present(tol))eps=tol
    allocate(model%centers(k,p),chosen(k));chosen=0
    do c=1,k
      do
        best=rng_integer(rng,1,n);if(.not.any(chosen(1:max(1,c-1))==best).or.c==1)exit
      end do
      chosen(c)=best;model%centers(c,:)=x(best,:)
    end do
    allocate(assign(n),countc(k),newc(k,p));model%converged=.false.
    do it=1,mi
      do i=1,n
        best=1;bestd=sum((x(i,:)-model%centers(1,:))**2)
        do c=2,k;d=sum((x(i,:)-model%centers(c,:))**2);if(d<bestd)then;best=c;bestd=d;end if;end do
        assign(i)=best
      end do
      newc=0.0_dp;countc=0
      do i=1,n;c=assign(i);newc(c,:)=newc(c,:)+x(i,:);countc(c)=countc(c)+1;end do
      do c=1,k
        if(countc(c)>0)then;newc(c,:)=newc(c,:)/real(countc(c),dp);else;newc(c,:)=model%centers(c,:);end if
      end do
      shift=maxval(abs(newc-model%centers));model%centers=newc
      if(shift<eps)then;model%converged=.true.;exit;end if
    end do
    model%iterations=min(it,mi);model%withinss=0.0_dp
    do i=1,n;model%withinss=model%withinss+sum((x(i,:)-model%centers(assign(i),:))**2);end do
  end subroutine
  subroutine predict_kmeans(model,x,response,prob)
    type(kmeans_model),intent(in)::model
    real(dp),intent(in)::x(:,:)
    integer,allocatable,intent(out)::response(:)
    real(dp),allocatable,intent(out),optional::prob(:,:)
    integer::i,c,best,k
    real(dp)::d,bestd,s
    real(dp),allocatable::w(:)
    k=size(model%centers,1);allocate(response(size(x,1)))
    if(present(prob))allocate(prob(size(x,1),k))
    allocate(w(k))
    do i=1,size(x,1)
      best=1;bestd=sum((x(i,:)-model%centers(1,:))**2)
      do c=2,k
        d=sum((x(i,:)-model%centers(c,:))**2)
        if(d<bestd)then;best=c;bestd=d;end if
      end do
      response(i)=best
      if(present(prob))then
        do c=1,k
          d=sqrt(sum((x(i,:)-model%centers(c,:))**2));w(c)=1.0_dp/max(d,1.0e-12_dp)
        end do
        s=sum(w);prob(i,:)=w/s
      end if
    end do
  end subroutine predict_kmeans

  subroutine fit_knn_regression(x,y,k,model)
    real(dp),intent(in)::x(:,:),y(:);integer,intent(in)::k;type(knn_model),intent(out)::model
    if(k<1.or.k>size(x,1))error stop 'fit_knn_regression: invalid k'
    model%x=x;model%y=y;model%k=k;model%classification=.false.
  end subroutine fit_knn_regression

  subroutine fit_knn_classifier(x,y,k,model)
    real(dp),intent(in)::x(:,:);integer,intent(in)::y(:),k;type(knn_model),intent(out)::model
    if(k<1.or.k>size(x,1))error stop 'fit_knn_classifier: invalid k'
    model%x=x;model%class_y=y;model%k=k;model%classification=.true.
  end subroutine fit_knn_classifier

  subroutine nearest_k(xtrain,q,k,idx)
    real(dp),intent(in)::xtrain(:,:),q(:);integer,intent(in)::k;integer,intent(out)::idx(k)
    real(dp),allocatable::d(:);integer::t,best
    allocate(d(size(xtrain,1)))
    do t=1,size(xtrain,1);d(t)=sum((xtrain(t,:)-q)**2);end do
    do t=1,k;best=minloc(d,dim=1);idx(t)=best;d(best)=huge(1.0_dp);end do
  end subroutine nearest_k

  subroutine predict_knn_regression(model,x,pred)
    type(knn_model),intent(in)::model;real(dp),intent(in)::x(:,:);real(dp),allocatable,intent(out)::pred(:)
    integer,allocatable::idx(:);integer::i
    if(model%classification)error stop 'predict_knn_regression: classifier model'
    allocate(pred(size(x,1)),idx(model%k))
    do i=1,size(x,1);call nearest_k(model%x,x(i,:),model%k,idx);pred(i)=sum(model%y(idx))/real(model%k,dp);end do
  end subroutine predict_knn_regression

  subroutine predict_knn_classifier(model,x,response,prob,nclass)
    type(knn_model),intent(in)::model;real(dp),intent(in)::x(:,:);integer,allocatable,intent(out)::response(:)
    real(dp),allocatable,intent(out),optional::prob(:,:);integer,intent(in),optional::nclass
    integer,allocatable::idx(:),counts(:);integer::i,j,c,nc
    if(.not.model%classification)error stop 'predict_knn_classifier: regression model'
    nc=maxval(model%class_y);if(present(nclass))nc=nclass
    allocate(response(size(x,1)),idx(model%k),counts(nc));if(present(prob))allocate(prob(size(x,1),nc))
    do i=1,size(x,1)
      call nearest_k(model%x,x(i,:),model%k,idx);counts=0
      do j=1,model%k;c=model%class_y(idx(j));if(c>=1.and.c<=nc)counts(c)=counts(c)+1;end do
      response(i)=maxloc(counts,dim=1)
      if(present(prob))prob(i,:)=real(counts,dp)/real(model%k,dp)
    end do
  end subroutine predict_knn_classifier
end module mlr_learners
