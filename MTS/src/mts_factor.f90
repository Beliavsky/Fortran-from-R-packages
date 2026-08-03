! SPDX-License-Identifier: Artistic-2.0
module mts_factor
   use mts_kinds, only : dp
   use mts_types, only : factor_result, bvar_result, mts_success, mts_invalid_input
   use mts_linalg, only : symmetric_eigen, matrix_sqrt_symmetric, inverse_matrix, least_squares, kronecker_product
   use mts_stats, only : column_mean, center_columns, covariance_matrix
   implicit none
   private

   public :: principal_components, asymptotic_pca, constrained_factor_model
   public :: stock_watson_forecast, fit_bvar

contains

   subroutine principal_components(x,n_factors,result,standardize)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::n_factors
      type(factor_result),intent(out)::result
      logical,intent(in),optional::standardize
      real(dp),allocatable::z(:,:),cov(:,:),vectors(:,:),values(:),mu(:),sd(:)
      logical::stand
      integer::n,k,r,j,istat
      n=size(x,1);k=size(x,2);r=min(max(1,n_factors),min(n,k));stand=.false.;if(present(standardize))stand=standardize
      if(n<2.or.k<1)then;result%status=mts_invalid_input;return;end if
      z=center_columns(x)
      if(stand)then
         allocate(sd(k));sd=sqrt(max(1.0e-14_dp,diagonal(covariance_matrix(x))))
         do j=1,k;z(:,j)=z(:,j)/sd(j);end do
      end if
      cov=matmul(transpose(z),z)/real(n,dp)
      call symmetric_eigen(cov,values,vectors,istat)
      result%loadings=vectors(:,1:r)
      result%factors=matmul(z,result%loadings)
      result%eigenvalues=values
      result%residuals=z-matmul(result%factors,transpose(result%loadings))
      result%status=istat
   end subroutine principal_components

   subroutine asymptotic_pca(x,n_factors,result)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::n_factors
      type(factor_result),intent(out)::result
      real(dp),allocatable::z(:,:),small_cov(:,:),vectors(:,:),values(:),load(:,:)
      integer::n,k,r,j,istat
      n=size(x,1);k=size(x,2);r=min(max(1,n_factors),min(n,k))
      if(n<2.or.k<1)then;result%status=mts_invalid_input;return;end if
      z=center_columns(x)
      if(k>n)then
         small_cov=matmul(z,transpose(z))/real(k,dp)
         call symmetric_eigen(small_cov,values,vectors,istat)
         result%factors=vectors(:,1:r)*sqrt(real(n,dp))
         allocate(load(k,r))
         do j=1,r
            if(values(j)>1.0e-14_dp)then
               load(:,j)=matmul(transpose(z),vectors(:,j))/sqrt(real(k,dp)*values(j))
            else
               load(:,j)=0.0_dp
            end if
         end do
         result%loadings=load
         result%eigenvalues=values
      else
         call principal_components(x,r,result)
         return
      end if
      result%residuals=z-matmul(result%factors,transpose(result%loadings))/sqrt(real(n,dp))
      result%status=istat
   end subroutine asymptotic_pca

   subroutine constrained_factor_model(x,h,n_factors,omega,factors,psi,loadings,status)
      real(dp),intent(in)::x(:,:),h(:,:)
      integer,intent(in)::n_factors
      real(dp),allocatable,intent(out)::omega(:,:),factors(:,:),psi(:,:),loadings(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::z(:,:),sd(:),hph(:,:),hpinvroot(:,:),yh(:,:),d(:,:),values(:),vectors(:,:),hp_inv(:,:),v1(:,:)
      integer::n,k,m,r,j,istat
      n=size(x,1);k=size(x,2);m=size(h,2);r=min(max(1,n_factors),m)
      if(size(h,1)/=k.or.n<2.or.m<1)then
         allocate(omega(0,0),factors(0,0),psi(0,0),loadings(0,0));if(present(status))status=mts_invalid_input;return
      end if
      z=center_columns(x);sd=sqrt(max(1.0e-14_dp,diagonal(covariance_matrix(x))))
      do j=1,k;z(:,j)=z(:,j)/sd(j);end do
      hph=matmul(transpose(h),h)
      call matrix_sqrt_symmetric(hph,hpinvroot,istat,inverse=.true.)
      yh=matmul(matmul(z,h),hpinvroot)
      d=matmul(transpose(yh),yh)/real(n,dp)
      call symmetric_eigen(d,values,vectors,istat)
      factors=matmul(yh,vectors(:,1:r))
      do j=1,r
         factors(:,j)=factors(:,j)/sqrt(max(1.0e-14_dp,sum((factors(:,j)-sum(factors(:,j))/real(n,dp))**2)/real(n-1,dp)))
      end do
      hp_inv=matmul(hpinvroot,hpinvroot)
      omega=matmul(hp_inv,matmul(transpose(matmul(z,h)),factors))/real(n,dp)
      loadings=matmul(h,omega)
      v1=covariance_matrix(z)
      psi=v1-matmul(loadings,transpose(loadings))
      if(present(status))status=istat
   end subroutine constrained_factor_model

   subroutine stock_watson_forecast(y,x,origin,n_factors,coefficients,forecast,mse,loadings,factors,status)
      real(dp),intent(in)::y(:),x(:,:)
      integer,intent(in)::origin,n_factors
      real(dp),allocatable,intent(out)::coefficients(:),forecast(:),loadings(:,:),factors(:,:)
      real(dp),intent(out)::mse
      integer,intent(out),optional::status
      type(factor_result)::pca
      real(dp),allocatable::z(:,:),mu(:),sd(:),design(:,:),beta(:,:),res(:,:)
      integer::n,k,r,j,istat,orig
      n=min(size(y),size(x,1));k=size(x,2);orig=min(max(2,origin),n);r=min(max(1,n_factors),k)
      allocate(mu(k),sd(k),z(n,k));mu=sum(x(1:orig,:),dim=1)/real(orig,dp)
      do j=1,k
         sd(j)=sqrt(max(1.0e-14_dp,sum((x(1:orig,j)-mu(j))**2)/real(orig-1,dp)))
         z(:,j)=(x(1:n,j)-mu(j))/sd(j)
      end do
      call principal_components(z(1:orig,:),r,pca)
      loadings=pca%loadings;factors=matmul(z,loadings)
      allocate(design(orig,r+1));design(:,1)=1.0_dp;design(:,2:)=factors(1:orig,:)
      call least_squares(design,reshape(y(1:orig),[orig,1]),beta,res,status=istat)
      coefficients=beta(:,1)
      if(orig<n)then
         forecast=coefficients(1)+matmul(factors(orig+1:n,:),coefficients(2:))
         mse=sum((y(orig+1:n)-forecast)**2)/real(n-orig,dp)
      else
         allocate(forecast(0));mse=0.0_dp
      end if
      if(present(status))status=istat
   end subroutine stock_watson_forecast

   subroutine fit_bvar(z,p,prior_precision,prior_scale,prior_df,result,prior_mean,include_mean)
      real(dp),intent(in)::z(:,:),prior_precision(:,:),prior_scale(:,:)
      integer,intent(in)::p,prior_df
      type(bvar_result),intent(out)::result
      real(dp),intent(in),optional::prior_mean(:,:)
      logical,intent(in),optional::include_mean
      real(dp),allocatable::x(:,:),y(:,:),phi0(:,:),w(:,:),winv(:,:),rhs(:,:),beta(:,:),res(:,:),diff(:,:),s(:,:)
      logical::mean_flag
      integer::n,k,d,lag,offset,istat,ne
      n=size(z,1);k=size(z,2);mean_flag=.true.;if(present(include_mean))mean_flag=include_mean
      ne=n-p;d=k*p+merge(1,0,mean_flag)
      if(p<1.or.ne<1.or.any(shape(prior_precision)/=[d,d]).or.any(shape(prior_scale)/=[k,k]))then
         result%status=mts_invalid_input;return
      end if
      allocate(x(ne,d),y(ne,k),phi0(d,k));x=0.0_dp;y=z(p+1:n,:);phi0=0.0_dp;offset=0
      if(mean_flag)then;x(:,1)=1.0_dp;offset=1;end if
      do lag=1,p;x(:,offset+(lag-1)*k+1:offset+lag*k)=z(p+1-lag:n-lag,:);end do
      if(present(prior_mean))then
         if(all(shape(prior_mean)==[d,k]))phi0=prior_mean
      end if
      w=matmul(transpose(x),x)+prior_precision
      call inverse_matrix(w,winv,istat)
      rhs=matmul(transpose(x),y)+matmul(prior_precision,phi0)
      beta=matmul(winv,rhs);res=y-matmul(x,beta);diff=beta-phi0
      s=matmul(transpose(res),res)+matmul(transpose(diff),matmul(prior_precision,diff))
      result%posterior_mean=beta;result%posterior_covariance=winv
      result%sigma=(prior_scale+s)/real(max(1,prior_df+ne-k-1),dp)
      result%order=p;result%include_mean=mean_flag;result%status=istat
   end subroutine fit_bvar

   pure function diagonal(a) result(d)
      real(dp),intent(in)::a(:,:)
      real(dp)::d(min(size(a,1),size(a,2)))
      integer::i
      do i=1,size(d);d(i)=a(i,i);end do
   end function diagonal

end module mts_factor
