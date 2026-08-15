! Correlated/heteroscedastic Gaussian GAMLSS using the supplied nlme backend.
! This is the matrix-first counterpart of combining a NO location model with
! nlme correlation and variance-function structures.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_correlation_v04
   use gamlss_kinds, only : dp
   use gamlss_fit, only : GAMLSS_NO
   use gamlss_types, only : gamlss_result_t, GAMLSS_METHOD_RS
   use nlme_gls, only : fit_gls
   use nlme_types, only : correlation_spec, variance_spec, gls_result, nlme_control
   implicit none
   private
   public :: correlated_no_result_t, fit_gamlss_no_gls

   type,public :: correlated_no_result_t
      type(gamlss_result_t) :: model
      type(gls_result) :: gls
      integer :: status=0
   end type correlated_no_result_t
contains

   subroutine fit_gamlss_no_gls(y,x_mu,result,correlation,variance,time,group,var_covariate, &
      var_group,coordinates,control)
      real(dp),intent(in)::y(:),x_mu(:,:)
      type(correlated_no_result_t),intent(out)::result
      type(correlation_spec),intent(in),optional::correlation
      type(variance_spec),intent(in),optional::variance
      real(dp),intent(in),optional::time(:),var_covariate(:),coordinates(:,:)
      integer,intent(in),optional::group(:),var_group(:)
      type(nlme_control),intent(in),optional::control
      integer::n,p,k
      real(dp),allocatable::sd(:)

      n=size(y);p=size(x_mu,2)
      if(n<=p.or.size(x_mu,1)/=n)then;result%status=1;return;end if
      call fit_gls(y,x_mu,result%gls,correlation=correlation,variance=variance,time=time,group=group, &
         var_covariate=var_covariate,var_group=var_group,coordinates=coordinates,control=control)
      result%status=result%gls%status
      if(result%status/=0)return

      result%model%family=GAMLSS_NO
      result%model%method=GAMLSS_METHOD_RS
      result%model%status=0
      result%model%converged=result%gls%converged
      result%model%iterations=result%gls%iterations
      result%model%mu%coefficients=result%gls%beta
      result%model%mu%covariance=result%gls%beta_cov
      result%model%mu%fitted=result%gls%fitted
      result%model%mu%eta=result%gls%fitted
      result%model%mu%edf=real(p,dp)
      result%model%residuals=result%gls%residuals
      allocate(sd(n));sd=0.0_dp
      if(allocated(result%gls%covariance))then
         do k=1,n;sd(k)=sqrt(max(0.0_dp,result%gls%covariance(k,k)));end do
      else
         sd=result%gls%sigma
      end if
      result%model%sigma%fitted=sd
      result%model%sigma%eta=log(max(sd,tiny(1.0_dp)))
      allocate(result%model%sigma%coefficients(1));result%model%sigma%coefficients=log(max(result%gls%sigma,tiny(1.0_dp)))
      allocate(result%model%sigma%covariance(1,1));result%model%sigma%covariance=0.0_dp
      result%model%global_deviance=-2.0_dp*result%gls%log_likelihood
      result%model%penalized_deviance=result%model%global_deviance
      result%model%aic=result%gls%aic;result%model%sbc=result%gls%bic
      result%model%df_fit=real(p+1,dp)
      if(allocated(result%gls%correlation_parameters))then
         result%model%df_fit=result%model%df_fit+real(size(result%gls%correlation_parameters),dp)
      end if
      if(allocated(result%gls%variance_parameters))then
         result%model%df_fit=result%model%df_fit+real(size(result%gls%variance_parameters),dp)
      end if
      result%model%df_residual=real(n,dp)-result%model%df_fit
      allocate(result%model%case_deviance(n))
      result%model%case_deviance=(result%gls%residuals/max(sd,sqrt(tiny(1.0_dp))))**2
   end subroutine fit_gamlss_no_gls
end module gamlss_correlation_v04
