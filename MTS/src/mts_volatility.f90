! SPDX-License-Identifier: Artistic-2.0
module mts_volatility
   use mts_kinds, only : dp
   use mts_types, only : dcc_model, bekk_model, mts_success, mts_invalid_input
   use mts_linalg, only : outer_product, cholesky_lower, nearest_psd, eye, symmetric_eigen, matrix_sqrt_symmetric
   use mts_stats, only : column_mean, center_columns, covariance_matrix, correlation_matrix, &
      multivariate_normal_logpdf, multivariate_student_t_logpdf
   use mts_optimize, only : bfgs_minimize
   use mts_regression, only : recursive_least_squares
   implicit none
   private

   public :: ewma_covariance, fit_ewma_lambda
   public :: dcc_correlations, dcc_log_likelihood, fit_dcc
   public :: bekk11_filter, bekk11_log_likelihood, fit_bekk11
   public :: modified_cholesky_volatility, common_volatility_components
   public :: constrained_group_correlation

contains

   subroutine ewma_covariance(returns,lambda,covariances,center,status)
      real(dp),intent(in)::returns(:,:),lambda
      real(dp),allocatable,intent(out)::covariances(:,:,:)
      logical,intent(in),optional::center
      integer,intent(out),optional::status
      real(dp),allocatable::x(:,:)
      logical::do_center
      integer::n,k,t
      n=size(returns,1);k=size(returns,2);do_center=.true.;if(present(center))do_center=center
      allocate(covariances(k,k,n))
      if(n<2.or.k<1.or.lambda<=0.0_dp.or.lambda>=1.0_dp)then
         covariances=0.0_dp;if(present(status))status=mts_invalid_input;return
      end if
      if(do_center)then;x=center_columns(returns);else;x=returns;end if
      covariances(:,:,1)=covariance_matrix(x)
      do t=2,n
         covariances(:,:,t)=lambda*covariances(:,:,t-1)+(1.0_dp-lambda)*outer_product(x(t-1,:),x(t-1,:))
         covariances(:,:,t)=0.5_dp*(covariances(:,:,t)+transpose(covariances(:,:,t)))
      end do
      if(present(status))status=mts_success
   end subroutine ewma_covariance

   subroutine fit_ewma_lambda(returns,lambda,covariances,log_likelihood,status,iterations)
      real(dp),intent(in)::returns(:,:)
      real(dp),intent(out)::lambda,log_likelihood
      real(dp),allocatable,intent(out)::covariances(:,:,:)
      integer,intent(out),optional::status,iterations
      real(dp)::par(1),lo(1),hi(1),value
      integer::istat,it
      par=0.96_dp;lo=1.0e-5_dp;hi=0.99999_dp
      call bfgs_minimize(objective,par,value,istat,it,200,1.0e-8_dp,lo,hi)
      lambda=par(1);log_likelihood=-value
      call ewma_covariance(returns,lambda,covariances)
      if(present(status))status=istat;if(present(iterations))iterations=it
   contains
      function objective(p) result(v)
         real(dp),intent(in)::p(:)
         real(dp)::v
         real(dp),allocatable::covs(:,:,:),x(:,:)
         integer::t,st
         x=center_columns(returns);call ewma_covariance(returns,p(1),covs)
         v=0.0_dp
         do t=2,size(returns,1)
            v=v-multivariate_normal_logpdf(x(t,:),0.0_dp*x(t,:),covs(:,:,t),st)
            if(st/=mts_success)then;v=huge(1.0_dp)/100.0_dp;return;end if
         end do
      end function objective
   end subroutine fit_ewma_lambda

   subroutine dcc_correlations(standardized_returns,alpha,beta,correlations,model_type,window,status)
      real(dp),intent(in)::standardized_returns(:,:),alpha,beta
      real(dp),allocatable,intent(out)::correlations(:,:,:)
      character(len=*),intent(in),optional::model_type
      integer,intent(in),optional::window
      integer,intent(out),optional::status
      real(dp),allocatable::qbar(:,:),q(:,:),local(:,:)
      real(dp)::d(size(standardized_returns,2))
      integer::n,k,t,m,i,j,istat
      character(len=16)::kind
      n=size(standardized_returns,1);k=size(standardized_returns,2);kind='engle';if(present(model_type))kind=adjustl(model_type)
      m=k+1;if(present(window))m=max(2,window)
      allocate(correlations(k,k,n));correlations=0.0_dp
      if(n<2.or.k<1.or.alpha<0.0_dp.or.beta<0.0_dp.or.alpha+beta>=1.0_dp)then
         if(present(status))status=mts_invalid_input;return
      end if
      qbar=correlation_matrix(standardized_returns);q=qbar;correlations(:,:,1)=qbar
      do t=2,n
         if(kind(1:1)=='t'.or.kind(1:1)=='T')then
            if(t>m)then
               local=correlation_matrix(standardized_returns(t-m:t-1,:))
            else
               local=qbar
            end if
            q=(1.0_dp-alpha-beta)*qbar+alpha*q+beta*local
         else
            q=(1.0_dp-alpha-beta)*qbar+alpha*outer_product(standardized_returns(t-1,:),standardized_returns(t-1,:))+beta*q
         end if
         d=sqrt(max(1.0e-14_dp,diagonal(q)))
         do i=1,k
            do j=1,k
               correlations(i,j,t)=q(i,j)/(d(i)*d(j))
            end do
            correlations(i,i,t)=1.0_dp
         end do
         call nearest_psd(correlations(:,:,t),local,1.0e-10_dp,istat)
         correlations(:,:,t)=correlation_from_covariance(local)
      end do
      if(present(status))status=mts_success
   end subroutine dcc_correlations

   function dcc_log_likelihood(standardized_returns,alpha,beta,model_type,distribution,degrees_freedom,window,correlations) result(loglik)
      real(dp),intent(in)::standardized_returns(:,:),alpha,beta
      character(len=*),intent(in),optional::model_type,distribution
      real(dp),intent(in),optional::degrees_freedom
      integer,intent(in),optional::window
      real(dp),allocatable,intent(out),optional::correlations(:,:,:)
      real(dp)::loglik,nu
      real(dp),allocatable::rho(:,:,:)
      character(len=16)::dist
      integer::t,istat
      dist='normal';if(present(distribution))dist=adjustl(distribution);nu=8.0_dp;if(present(degrees_freedom))nu=degrees_freedom
      call dcc_correlations(standardized_returns,alpha,beta,rho,model_type,window,istat)
      if(istat/=mts_success)then;loglik=-huge(1.0_dp);return;end if
      loglik=0.0_dp
      do t=2,size(standardized_returns,1)
         if(dist(1:1)=='t'.or.dist(1:1)=='T'.or.dist(1:1)=='s'.or.dist(1:1)=='S')then
            loglik=loglik+multivariate_student_t_logpdf(standardized_returns(t,:),0.0_dp*standardized_returns(t,:),rho(:,:,t),nu,.true.,istat)
         else
            loglik=loglik+multivariate_normal_logpdf(standardized_returns(t,:),0.0_dp*standardized_returns(t,:),rho(:,:,t),istat)
         end if
         if(istat/=mts_success)then;loglik=-huge(1.0_dp);exit;end if
      end do
      if(present(correlations))correlations=rho
   end function dcc_log_likelihood

   subroutine fit_dcc(standardized_returns,model,model_type,distribution,initial,window,max_iterations,tolerance)
      real(dp),intent(in)::standardized_returns(:,:)
      type(dcc_model),intent(out)::model
      character(len=*),intent(in),optional::model_type,distribution
      real(dp),intent(in),optional::initial(:)
      integer,intent(in),optional::window,max_iterations
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::par(:),lo(:),hi(:)
      real(dp)::value,tol
      integer::np,it,istat,maxit,m
      character(len=16)::kind,dist
      kind='engle';if(present(model_type))kind=adjustl(model_type)
      dist='normal';if(present(distribution))dist=adjustl(distribution)
      np=2;if(dist(1:1)=='t'.or.dist(1:1)=='T'.or.dist(1:1)=='s'.or.dist(1:1)=='S')np=3
      allocate(par(np),lo(np),hi(np));par(1:2)=[0.02_dp,0.95_dp]
      if(present(initial))par(1:min(np,size(initial)))=initial(1:min(np,size(initial)))
      lo=1.0e-6_dp;hi=0.999_dp
      if(np==3)then
         par(3)=8.0_dp
         if(present(initial))then
            if(size(initial)>=3)par(3)=initial(3)
         end if
         par(3)=max(4.1_dp,par(3));lo(3)=4.01_dp;hi(3)=100.0_dp
      end if
      maxit=300;if(present(max_iterations))maxit=max_iterations;tol=1.0e-7_dp;if(present(tolerance))tol=tolerance
      m=size(standardized_returns,2)+1;if(present(window))m=window
      call bfgs_minimize(objective,par,value,istat,it,maxit,tol,lo,hi)
      model%model_type=kind;model%distribution=dist;model%alpha=par(1);model%beta=par(2)
      if(np==3)model%degrees_freedom=par(3)
      model%log_likelihood=-value;model%iterations=it;model%status=istat
      model%unconditional_corr=correlation_matrix(standardized_returns)
      call dcc_correlations(standardized_returns,model%alpha,model%beta,model%correlations,kind,m)
   contains
      function objective(p) result(v)
         real(dp),intent(in)::p(:)
         real(dp)::v,ll
         if(p(1)+p(2)>=0.9999_dp)then;v=1.0e8_dp+1.0e8_dp*(p(1)+p(2)-0.9999_dp)**2;return;end if
         if(size(p)==3)then
            ll=dcc_log_likelihood(standardized_returns,p(1),p(2),kind,dist,p(3),m)
         else
            ll=dcc_log_likelihood(standardized_returns,p(1),p(2),kind,dist,window=m)
         end if
         if(ll<=-huge(1.0_dp)/2.0_dp)then;v=huge(1.0_dp)/100.0_dp;else;v=-ll;end if
      end function objective
   end subroutine fit_dcc

   subroutine bekk11_filter(returns,mean,c,a,b,covariances,residuals,status)
      real(dp),intent(in)::returns(:,:),mean(:),c(:,:),a(:,:),b(:,:)
      real(dp),allocatable,intent(out)::covariances(:,:,:),residuals(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::tmp(:,:)
      integer::n,k,t,istat,overall
      n=size(returns,1);k=size(returns,2);allocate(covariances(k,k,n),residuals(n,k));overall=mts_success
      if(size(mean)/=k.or.any(shape(c)/=[k,k]).or.any(shape(a)/=[k,k]).or.any(shape(b)/=[k,k]))then
         covariances=0.0_dp;residuals=0.0_dp;if(present(status))status=mts_invalid_input;return
      end if
      residuals=returns-spread(mean,dim=1,ncopies=n);covariances(:,:,1)=covariance_matrix(residuals)
      do t=2,n
         covariances(:,:,t)=matmul(c,transpose(c))+matmul(a,matmul(outer_product(residuals(t-1,:),residuals(t-1,:)),transpose(a)))+ &
            matmul(b,matmul(covariances(:,:,t-1),transpose(b)))
         covariances(:,:,t)=0.5_dp*(covariances(:,:,t)+transpose(covariances(:,:,t)))
         call nearest_psd(covariances(:,:,t),tmp,1.0e-12_dp,istat);covariances(:,:,t)=tmp
         if(istat/=mts_success)overall=istat
      end do
      if(present(status))status=overall
   end subroutine bekk11_filter

   function bekk11_log_likelihood(returns,mean,c,a,b,covariances) result(loglik)
      real(dp),intent(in)::returns(:,:),mean(:),c(:,:),a(:,:),b(:,:)
      real(dp),allocatable,intent(out),optional::covariances(:,:,:)
      real(dp)::loglik
      real(dp),allocatable::covs(:,:,:),res(:,:)
      integer::t,istat
      call bekk11_filter(returns,mean,c,a,b,covs,res,istat)
      if(istat/=mts_success)then;loglik=-huge(1.0_dp);return;end if
      loglik=0.0_dp
      do t=2,size(returns,1)
         loglik=loglik+multivariate_normal_logpdf(res(t,:),0.0_dp*res(t,:),covs(:,:,t),istat)
         if(istat/=mts_success)then;loglik=-huge(1.0_dp);exit;end if
      end do
      if(present(covariances))covariances=covs
   end function bekk11_log_likelihood

   subroutine fit_bekk11(returns,model,include_mean,max_iterations,tolerance)
      real(dp),intent(in)::returns(:,:)
      type(bekk_model),intent(out)::model
      logical,intent(in),optional::include_mean
      integer,intent(in),optional::max_iterations
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::par(:),lo(:),hi(:),l(:,:),cov(:,:),mean(:),c(:,:),a(:,:),b(:,:),res(:,:)
      real(dp)::value,tol
      logical::mean_flag
      integer::n,k,nc,np,offset,i,j,it,istat,maxit
      n=size(returns,1);k=size(returns,2);mean_flag=.true.;if(present(include_mean))mean_flag=include_mean
      nc=k*(k+1)/2;np=merge(k,0,mean_flag)+nc+2*k*k
      allocate(par(np),lo(np),hi(np));par=0.0_dp;lo=-0.999_dp;hi=0.999_dp;offset=0
      if(mean_flag)then;par(1:k)=column_mean(returns);lo(1:k)=minval(returns,dim=1);hi(1:k)=maxval(returns,dim=1);offset=k;end if
      cov=covariance_matrix(returns);call cholesky_lower(cov,l,istat,jitter=1.0e-10_dp)
      do j=1,k
         do i=j,k
            offset=offset+1;par(offset)=0.25_dp*l(i,j)
            if(i==j)then;lo(offset)=1.0e-8_dp;hi(offset)=max(1.0_dp,2.0_dp*l(i,j));end if
         end do
      end do
      do i=1,k;par(offset+(i-1)*k+i)=0.10_dp;end do;offset=offset+k*k
      do i=1,k;par(offset+(i-1)*k+i)=0.85_dp;end do
      maxit=400;if(present(max_iterations))maxit=max_iterations;tol=1.0e-6_dp;if(present(tolerance))tol=tolerance
      call bfgs_minimize(objective,par,value,istat,it,maxit,tol,lo,hi)
      call unpack_bekk(par,k,mean_flag,mean,c,a,b)
      model%mean=mean;model%c=c;model%a=a;model%b=b;model%include_mean=mean_flag
      model%log_likelihood=-value;model%iterations=it;model%status=istat
      call bekk11_filter(returns,mean,c,a,b,model%covariance,res)
   contains
      function objective(p) result(v)
         real(dp),intent(in)::p(:)
         real(dp)::v,ll,pen
         real(dp),allocatable::m0(:),c0(:,:),a0(:,:),b0(:,:)
         call unpack_bekk(p,k,mean_flag,m0,c0,a0,b0)
         pen=max(0.0_dp,frobenius(a0)+frobenius(b0)-1.8_dp)
         ll=bekk11_log_likelihood(returns,m0,c0,a0,b0)
         if(ll<=-huge(1.0_dp)/2.0_dp)then;v=huge(1.0_dp)/100.0_dp;else;v=-ll+1.0e4_dp*pen*pen;end if
      end function objective
   end subroutine fit_bekk11

   subroutine unpack_bekk(par,k,include_mean,mean,c,a,b)
      real(dp),intent(in)::par(:)
      integer,intent(in)::k
      logical,intent(in)::include_mean
      real(dp),allocatable,intent(out)::mean(:),c(:,:),a(:,:),b(:,:)
      integer::offset,i,j
      allocate(mean(k),c(k,k),a(k,k),b(k,k));mean=0.0_dp;c=0.0_dp;offset=0
      if(include_mean)then;mean=par(1:k);offset=k;end if
      do j=1,k
         do i=j,k;offset=offset+1;c(i,j)=par(offset);end do
      end do
      a=reshape(par(offset+1:offset+k*k),[k,k]);offset=offset+k*k
      b=reshape(par(offset+1:offset+k*k),[k,k])
   end subroutine unpack_bekk

   subroutine modified_cholesky_volatility(returns,window,lambda,covariances,coefficients,components,status)
      real(dp),intent(in)::returns(:,:),lambda
      integer,intent(in)::window
      real(dp),allocatable,intent(out)::covariances(:,:,:),coefficients(:,:),components(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::centered(:,:),path(:,:),coef(:),res(:),a(:,:),ainv(:,:),vol(:,:)
      integer::n,k,i,t,count,istat
      n=size(returns,1);k=size(returns,2);centered=center_columns(returns)
      allocate(coefficients(max(0,n-window),k*(k-1)/2),components(max(0,n-window),k),vol(max(0,n-window),k))
      coefficients=0.0_dp;components=0.0_dp;count=0
      components(:,1)=centered(window+1:n,1)
      do i=2,k
         call recursive_least_squares(centered(:,i),centered(:,1:i-1),window,coef,path,lambda,istat)
         coefficients(:,count+1:count+i-1)=-path(window+1:n,:)
         do t=window+1,n
            components(t-window,i)=centered(t,i)-dot_product(path(t,:),centered(t,1:i-1))
         end do
         count=count+i-1
      end do
      do i=1,k
         vol(1,i)=max(1.0e-10_dp,sum(components(:,i)**2)/real(max(1,size(components,1)),dp))
         do t=2,size(components,1);vol(t,i)=lambda*vol(t-1,i)+(1.0_dp-lambda)*components(t-1,i)**2;end do
      end do
      allocate(covariances(k,k,size(components,1)))
      do t=1,size(components,1)
         a=eye(k);count=0
         do i=2,k;a(i,1:i-1)=coefficients(t,count+1:count+i-1);count=count+i-1;end do
         call invert_lower_unit(a,ainv)
         covariances(:,:,t)=matmul(ainv,matmul(diagonal_matrix(vol(t,:)),transpose(ainv)))
      end do
      if(present(status))status=mts_success
   end subroutine modified_cholesky_volatility

   subroutine common_volatility_components(returns,max_lag,loadings,eigenvalues,transformed,standardize,status)
      real(dp),intent(in)::returns(:,:)
      integer,intent(in)::max_lag
      real(dp),allocatable,intent(out)::loadings(:,:),eigenvalues(:),transformed(:,:)
      logical,intent(in),optional::standardize
      integer,intent(out),optional::status
      real(dp),allocatable::x(:,:),root_inv(:,:),a(:,:),cmtx(:,:),vectors(:,:)
      real(dp)::y1(size(returns,1)-max_lag),y2(size(returns,1)-max_lag)
      logical::stand
      integer::k,n,h,i,j,ii,jj,istat
      n=size(returns,1);k=size(returns,2);stand=.false.;if(present(standardize))stand=standardize
      if(max_lag<1.or.n<=max_lag)then;if(present(status))status=mts_invalid_input;return;end if
      x=center_columns(returns);call matrix_sqrt_symmetric(covariance_matrix(x),root_inv,istat,inverse=.true.)
      x=matmul(x,root_inv);allocate(a(k,k),cmtx(k,k));a=0.0_dp
      do h=1,max_lag
         do i=1,k
            do j=i,k
               y2=x(max_lag+1-h:n-h,i)*x(max_lag+1-h:n-h,j);cmtx=0.0_dp
               do ii=1,k
                  do jj=ii,k
                     y1=x(max_lag+1:n,ii)*x(max_lag+1:n,jj)
                     cmtx(ii,jj)=cov_pair(y1,y2);cmtx(jj,ii)=cmtx(ii,jj)
                  end do
               end do
               a=a+matmul(cmtx,cmtx)
            end do
         end do
      end do
      if(stand)a=correlation_from_covariance(a)
      call symmetric_eigen(a,eigenvalues,vectors,istat);loadings=matmul(root_inv,vectors);transformed=matmul(returns,loadings)
      if(present(status))status=istat
   end subroutine common_volatility_components

   subroutine constrained_group_correlation(returns,group_sizes,end_index,span,unconstrained,constrained,status)
      !! Sample correlation with common within-group and between-group entries.
      real(dp),intent(in)::returns(:,:)
      integer,intent(in)::group_sizes(:)
      integer,intent(in),optional::end_index,span
      real(dp),allocatable,intent(out)::unconstrained(:,:),constrained(:,:)
      integer,intent(out),optional::status
      integer::n,k,last,width,first,ng,g1,g2,i,j,s1,e1,s2,e2,count
      real(dp)::average
      n=size(returns,1);k=size(returns,2);ng=size(group_sizes)
      last=n;if(present(end_index))last=min(n,max(1,end_index))
      width=last;if(present(span))width=max(2,min(last,span));first=last-width+1
      if(ng<1.or.any(group_sizes<1).or.sum(group_sizes)/=k.or.width<2)then
         allocate(unconstrained(0,0),constrained(0,0));if(present(status))status=mts_invalid_input;return
      end if
      unconstrained=correlation_matrix(returns(first:last,:));constrained=unconstrained
      s1=1
      do g1=1,ng
         e1=s1+group_sizes(g1)-1
         if(group_sizes(g1)>1)then
            average=0.0_dp;count=0
            do i=s1,e1-1
               do j=i+1,e1;average=average+unconstrained(i,j);count=count+1;end do
            end do
            average=average/real(count,dp)
            do i=s1,e1-1
               do j=i+1,e1;constrained(i,j)=average;constrained(j,i)=average;end do
            end do
         end if
         s2=e1+1
         do g2=g1+1,ng
            e2=s2+group_sizes(g2)-1;average=0.0_dp;count=0
            do i=s1,e1
               do j=s2,e2;average=average+unconstrained(i,j);count=count+1;end do
            end do
            average=average/real(count,dp)
            do i=s1,e1
               do j=s2,e2;constrained(i,j)=average;constrained(j,i)=average;end do
            end do
            s2=e2+1
         end do
         s1=e1+1
      end do
      if(present(status))status=mts_success
   end subroutine constrained_group_correlation

   function correlation_from_covariance(cov) result(corr)
      real(dp),intent(in)::cov(:,:)
      real(dp)::corr(size(cov,1),size(cov,2)),sd(size(cov,1))
      integer::i,j
      sd=sqrt(max(1.0e-14_dp,diagonal(cov)))
      do i=1,size(cov,1)
         do j=1,size(cov,2);corr(i,j)=cov(i,j)/(sd(i)*sd(j));end do
         corr(i,i)=1.0_dp
      end do
   end function correlation_from_covariance

   pure function diagonal(a) result(d)
      real(dp),intent(in)::a(:,:)
      real(dp)::d(min(size(a,1),size(a,2)))
      integer::i
      do i=1,size(d);d(i)=a(i,i);end do
   end function diagonal

   pure function diagonal_matrix(d) result(a)
      real(dp),intent(in)::d(:)
      real(dp)::a(size(d),size(d))
      integer::i
      a=0.0_dp;do i=1,size(d);a(i,i)=d(i);end do
   end function diagonal_matrix

   pure function frobenius(a) result(v)
      real(dp),intent(in)::a(:,:)
      real(dp)::v
      v=sqrt(sum(a*a))
   end function frobenius

   subroutine invert_lower_unit(a,ainv)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable,intent(out)::ainv(:,:)
      integer::n,i,j,k
      n=size(a,1);allocate(ainv(n,n));ainv=0.0_dp
      do i=1,n;ainv(i,i)=1.0_dp;end do
      do i=2,n
         do j=1,i-1
            do k=j,i-1;ainv(i,j)=ainv(i,j)-a(i,k)*ainv(k,j);end do
         end do
      end do
   end subroutine invert_lower_unit

   pure function cov_pair(x,y) result(c)
      real(dp),intent(in)::x(:),y(:)
      real(dp)::c,mx,my
      mx=sum(x)/real(size(x),dp);my=sum(y)/real(size(y),dp)
      c=sum((x-mx)*(y-my))/real(max(1,size(x)-1),dp)
   end function cov_pair

end module mts_volatility
