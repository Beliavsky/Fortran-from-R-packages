! SPDX-License-Identifier: Artistic-2.0
module mts_diagnostics
   use mts_kinds, only : dp
   use mts_types, only : diagnostic_result, mts_success, mts_invalid_input
   use mts_linalg, only : inverse_matrix, kronecker_product, vec_matrix, vech_matrix, least_squares
   use mts_stats, only : covariance_matrix, autocovariance_matrix, correlation_matrix, chi_square_survival, ranks
   implicit none
   private

   public :: cross_correlation_matrices, multivariate_portmanteau
   public :: multivariate_arch_test, rank_arch_test, volatility_diagnostics
   public :: standardized_residuals

contains

   subroutine cross_correlation_matrices(x,max_lag,ccm,p_values,status)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::max_lag
      real(dp),allocatable,intent(out)::ccm(:,:,:),p_values(:)
      integer,intent(out),optional::status
      real(dp),allocatable::corr0_inv(:,:),kron_inv(:,:),v(:)
      real(dp)::stat
      integer::k,n,h,istat
      n=size(x,1);k=size(x,2)
      if(max_lag<0.or.n<=max_lag.or.k<1)then
         allocate(ccm(k,k,0),p_values(0));if(present(status))status=mts_invalid_input;return
      end if
      allocate(ccm(k,k,0:max_lag),p_values(max_lag));ccm(:,:,0)=correlation_matrix(x)
      call inverse_matrix(ccm(:,:,0),corr0_inv,istat)
      if(istat/=mts_success)then;if(present(status))status=istat;return;end if
      kron_inv=kronecker_product(corr0_inv,corr0_inv)
      do h=1,max_lag
         ccm(:,:,h)=cross_corr_lag(x,h)
         v=vec_matrix(ccm(:,:,h))
         stat=real(n*n,dp)/real(n-h,dp)*dot_product(v,matmul(kron_inv,v))
         p_values(h)=chi_square_survival(max(0.0_dp,stat),k*k)
      end do
      if(present(status))status=mts_success
   end subroutine cross_correlation_matrices

   function cross_corr_lag(x,lag) result(corr)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::lag
      real(dp)::corr(size(x,2),size(x,2)),mu(size(x,2)),sd(size(x,2))
      real(dp),allocatable::y(:,:)
      integer::n,i,j
      n=size(x,1);mu=sum(x,dim=1)/real(n,dp);allocate(y(n,size(x,2)))
      do i=1,n;y(i,:)=x(i,:)-mu;end do
      sd=sqrt(max(1.0e-30_dp,sum(y*y,dim=1)/real(max(1,n-1),dp)))
      corr=0.0_dp
      do i=1,size(x,2)
         do j=1,size(x,2)
            corr(i,j)=dot_product(y(lag+1:n,i),y(1:n-lag,j))/real(n,dp)/(sd(i)*sd(j))
         end do
      end do
   end function cross_corr_lag

   subroutine multivariate_portmanteau(x,max_lag,result,model_order,adjusted)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::max_lag
      type(diagnostic_result),intent(out)::result
      integer,intent(in),optional::model_order
      logical,intent(in),optional::adjusted
      real(dp),allocatable::g0inv(:,:),gh(:,:)
      real(dp)::q,weight
      integer::n,k,h,adj,istat
      logical::use_adjusted
      n=size(x,1);k=size(x,2);adj=0;if(present(model_order))adj=max(0,model_order)
      use_adjusted=.false.;if(present(adjusted))use_adjusted=adjusted
      if(max_lag<1.or.n<=max_lag)then;result%status=mts_invalid_input;return;end if
      call inverse_matrix(autocovariance_matrix(x,0),g0inv,istat)
      if(istat/=mts_success)then;result%status=istat;return;end if
      q=0.0_dp
      do h=1,max_lag
         gh=autocovariance_matrix(x,h)
         weight=1.0_dp
         if(use_adjusted)weight=real(n+2,dp)/real(n-h,dp)
         q=q+weight*trace_local(matmul(transpose(gh),matmul(g0inv,matmul(gh,g0inv))))
      end do
      if(use_adjusted)then;q=real(n,dp)*q;else;q=real(n*n,dp)*q;end if
      result%statistic=max(0.0_dp,q)
      result%degrees_freedom=max(1,k*k*(max_lag-adj))
      result%p_value=chi_square_survival(result%statistic,result%degrees_freedom)
      result%status=mts_success
   end subroutine multivariate_portmanteau

   subroutine multivariate_arch_test(x,max_lag,result)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::max_lag
      type(diagnostic_result),intent(out)::result
      real(dp),allocatable::u(:,:),design(:,:),y(:,:),beta(:,:),res(:,:),covb(:,:)
      real(dp)::sst,sse,r2
      integer::n,k,m,t,j,d,istat
      n=size(x,1);k=size(x,2);m=k*(k+1)/2
      if(max_lag<1.or.n<=max_lag+1)then;result%status=mts_invalid_input;return;end if
      allocate(u(n,m))
      do t=1,n;u(t,:)=vech_matrix(outer_local(x(t,:),x(t,:)));end do
      allocate(design(n-max_lag,1+m*max_lag),y(n-max_lag,m));design(:,1)=1.0_dp;y=u(max_lag+1:n,:)
      do j=1,max_lag
         design(:,1+(j-1)*m+1:1+j*m)=u(max_lag+1-j:n-j,:)
      end do
      call least_squares(design,y,beta,res,covb,istat,ridge=1.0e-10_dp)
      sst=sum((y-spread(sum(y,dim=1)/real(size(y,1),dp),dim=1,ncopies=size(y,1)))**2)
      sse=sum(res**2);r2=max(0.0_dp,min(1.0_dp,1.0_dp-sse/max(sst,tiny(1.0_dp))))
      result%statistic=real(size(y,1),dp)*real(m,dp)*r2
      result%degrees_freedom=m*m*max_lag
      result%p_value=chi_square_survival(result%statistic,result%degrees_freedom)
      result%status=istat
   end subroutine multivariate_arch_test

   subroutine rank_arch_test(x,max_lag,result)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::max_lag
      type(diagnostic_result),intent(out)::result
      real(dp),allocatable::ranked(:,:),squares(:,:)
      integer::j,n
      n=size(x,1);allocate(ranked(n,size(x,2)),squares(n,size(x,2)))
      do j=1,size(x,2)
         ranked(:,j)=(ranks(x(:,j))-0.5_dp)/real(n,dp)-0.5_dp
      end do
      squares=ranked*ranked
      call multivariate_portmanteau(squares,max_lag,result,adjusted=.true.)
   end subroutine rank_arch_test

   subroutine standardized_residuals(residuals,covariances,z,status)
      use mts_linalg, only : matrix_sqrt_symmetric
      real(dp),intent(in)::residuals(:,:),covariances(:,:,:)
      real(dp),allocatable,intent(out)::z(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::root_inv(:,:)
      integer::t,n,k,istat,overall
      n=size(residuals,1);k=size(residuals,2);allocate(z(n,k));z=0.0_dp;overall=mts_success
      if(size(covariances,1)/=k.or.size(covariances,2)/=k.or.size(covariances,3)/=n)then
         if(present(status))status=mts_invalid_input;return
      end if
      do t=1,n
         call matrix_sqrt_symmetric(covariances(:,:,t),root_inv,istat,inverse=.true.)
         if(istat==mts_success)then;z(t,:)=matmul(root_inv,residuals(t,:));else;overall=istat;end if
      end do
      if(present(status))status=overall
   end subroutine standardized_residuals

   subroutine volatility_diagnostics(residuals,covariances,max_lag,level_test,squared_test,status)
      real(dp),intent(in)::residuals(:,:),covariances(:,:,:)
      integer,intent(in)::max_lag
      type(diagnostic_result),intent(out)::level_test,squared_test
      integer,intent(out),optional::status
      real(dp),allocatable::z(:,:)
      integer::istat
      call standardized_residuals(residuals,covariances,z,istat)
      call multivariate_portmanteau(z,max_lag,level_test)
      call multivariate_arch_test(z,max_lag,squared_test)
      if(present(status))status=istat
   end subroutine volatility_diagnostics

   pure function trace_local(a) result(v)
      real(dp),intent(in)::a(:,:)
      real(dp)::v
      integer::i
      v=0.0_dp
      do i=1,min(size(a,1),size(a,2));v=v+a(i,i);end do
   end function trace_local

   pure function outer_local(x,y) result(a)
      real(dp),intent(in)::x(:),y(:)
      real(dp)::a(size(x),size(y))
      integer::i
      do i=1,size(x);a(i,:)=x(i)*y;end do
   end function outer_local

end module mts_diagnostics
