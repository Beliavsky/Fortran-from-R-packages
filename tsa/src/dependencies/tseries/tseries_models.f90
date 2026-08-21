! Experimental modern Fortran translation of computational routines from
! the R package tseries 0.10-62. Original authors include Adrian Trapletti
! and Kurt Hornik; Blake LeBaron contributed the original BDS code.
! Licensed under GPL-2.0-only OR GPL-3.0-only. See LICENSE and NOTICE.

module tseries_models
   use tseries_kinds, only : dp
   use tseries_types, only : arma_result, garch_result
   use tseries_stats, only : mean_value, variance_value
   use tseries_optimize, only : nelder_mead, numerical_hessian
   use tseries_linalg, only : invert_matrix
   implicit none
   private

   public :: arma_fit
   public :: arma_residuals
   public :: garch_fit
   public :: garch_variance

contains

   subroutine arma_residuals(x, ar_lags, ma_lags, coefficients, include_intercept, residuals, fitted)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: ar_lags(:), ma_lags(:)
      real(dp), intent(in) :: coefficients(:)
      logical, intent(in) :: include_intercept
      real(dp), intent(out) :: residuals(:)
      real(dp), intent(out), optional :: fitted(:)
      integer :: i,j,na,nm,max_lag,intercept_index
      real(dp) :: prediction

      na=size(ar_lags); nm=size(ma_lags)
      max_lag=0
      if(na>0) max_lag=max(max_lag,maxval(ar_lags))
      if(nm>0) max_lag=max(max_lag,maxval(ma_lags))
      residuals=0.0_dp
      intercept_index=na+nm+1
      do i=max_lag+1,size(x)
         prediction=0.0_dp
         if(include_intercept) prediction=coefficients(intercept_index)
         do j=1,na
            prediction=prediction+coefficients(j)*x(i-ar_lags(j))
         end do
         do j=1,nm
            prediction=prediction+coefficients(na+j)*residuals(i-ma_lags(j))
         end do
         residuals(i)=x(i)-prediction
      end do
      if(present(fitted)) fitted=x-residuals
   end subroutine arma_residuals

   function arma_fit(x,p,q,include_intercept,max_iterations,tolerance) result(fit)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: p,q
      logical, intent(in), optional :: include_intercept
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(arma_result) :: fit
      integer, allocatable :: ar_lags(:),ma_lags(:)
      real(dp), allocatable :: par(:),hess(:,:),hinv(:,:),resid(:),fitted(:)
      integer :: ncoef,maxit,it,status,i,max_lag
      logical :: intercept
      real(dp) :: tol,obj,sigma2

      intercept=.true.; maxit=2000; tol=1.0e-8_dp
      if(present(include_intercept)) intercept=include_intercept
      if(present(max_iterations)) maxit=max_iterations
      if(present(tolerance)) tol=tolerance
      fit%p=p; fit%q=q; fit%include_intercept=intercept
      if(size(x)<=max(p,q)+2 .or. p<0 .or. q<0) then
         fit%status=1; fit%message='invalid ARMA dimensions'; return
      end if
      allocate(ar_lags(p),ma_lags(q))
      ar_lags=[(i,i=1,p)]; ma_lags=[(i,i=1,q)]
      ncoef=p+q+merge(1,0,intercept)
      allocate(par(ncoef),resid(size(x)),fitted(size(x)),hess(ncoef,ncoef),hinv(ncoef,ncoef))
      par=0.0_dp
      if(intercept) par(ncoef)=mean_value(x)
      call nelder_mead(objective,par,obj,it,status,max_iterations=maxit,tolerance=tol,step=0.05_dp)
      call arma_residuals(x,ar_lags,ma_lags,par,intercept,resid,fitted)
      max_lag=max(p,q)
      allocate(fit%coefficients(ncoef),fit%residuals(size(x)),fit%fitted(size(x)),fit%covariance(ncoef,ncoef))
      fit%coefficients=par; fit%residuals=resid; fit%fitted=fitted
      fit%css=sum(resid(max_lag+1:)**2)
      sigma2=fit%css/real(size(x)-max_lag,dp)
      fit%aic=real(size(x),dp)*(1.0_dp+log(2.0_dp*acos(-1.0_dp))+log(max(sigma2,tiny(1.0_dp))))+2.0_dp*real(ncoef,dp)
      fit%iterations=it; fit%status=status
      if(status==0) then
         fit%message='converged'
      else
         fit%message='optimizer reached iteration limit'
      end if
      call numerical_hessian(objective,par,hess,status)
      if(status==0) then
         call invert_matrix(hess,hinv,status)
      end if
      if(status==0) then
         fit%covariance=2.0_dp*fit%css/real(size(x),dp)*hinv
      else
         fit%covariance=0.0_dp
      end if

   contains
      function objective(a) result(value)
         real(dp), intent(in) :: a(:)
         real(dp) :: value
         real(dp) :: u(size(x))
         call arma_residuals(x,ar_lags,ma_lags,a,intercept,u)
         value=sum(u(max(p,q)+1:)**2)
         if(.not.(value<huge(value))) value=huge(value)/100.0_dp
      end function objective
   end function arma_fit

   subroutine garch_variance(x,coefficients,p,q,h,genuine)
      real(dp), intent(in) :: x(:),coefficients(:)
      integer, intent(in) :: p,q
      real(dp), intent(out) :: h(:)
      logical, intent(in), optional :: genuine
      logical :: predict_next
      integer :: nout,i,j,maxpq
      real(dp) :: unconditional,coef_sum,value

      predict_next=.false.
      if(present(genuine)) predict_next=genuine
      nout=size(x)+merge(1,0,predict_next)
      if(size(h)/=nout) return
      maxpq=max(p,q)
      coef_sum=sum(coefficients(2:))
      if(coef_sum<1.0_dp) then
         unconditional=coefficients(1)/max(1.0e-12_dp,1.0_dp-coef_sum)
      else
         unconditional=variance_value(x,.false.)
      end if
      h=max(unconditional,1.0e-12_dp)
      do i=maxpq+1,nout
         value=coefficients(1)
         do j=1,q
            value=value+coefficients(1+j)*x(i-j)**2
         end do
         do j=1,p
            value=value+coefficients(1+q+j)*h(i-j)
         end do
         h(i)=max(value,1.0e-12_dp)
      end do
   end subroutine garch_variance

   function garch_fit(x,p,q,max_iterations,tolerance) result(fit)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: p,q
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(garch_result) :: fit
      real(dp), allocatable :: theta(:),coef(:),h(:),resid(:),hess(:,:),hinv(:,:)
      integer :: ncoef,maxit,it,status,maxpq
      real(dp) :: tol,obj,v,small

      maxit=2000; tol=1.0e-8_dp
      if(present(max_iterations)) maxit=max_iterations
      if(present(tolerance)) tol=tolerance
      fit%p=p; fit%q=q
      if(size(x)<=max(p,q)+2 .or. p<0 .or. q<0) then
         fit%status=1; fit%message='invalid GARCH dimensions'; return
      end if
      ncoef=1+p+q; maxpq=max(p,q)
      allocate(theta(ncoef),coef(ncoef),h(size(x)),resid(size(x)),hess(ncoef,ncoef),hinv(ncoef,ncoef))
      v=max(variance_value(x,.false.),1.0e-8_dp)
      small=0.05_dp
      coef(1)=v*(1.0_dp-small*real(ncoef-1,dp))
      if(ncoef>1) coef(2:)=small
      call coefficients_to_theta(coef,theta)
      call nelder_mead(objective,theta,obj,it,status,max_iterations=maxit,tolerance=tol,step=0.10_dp)
      call theta_to_coefficients(theta,coef)
      call garch_variance(x,coef,p,q,h)
      resid=x/sqrt(h)
      allocate(fit%coefficients(ncoef),fit%conditional_variance(size(x)),fit%residuals(size(x)),fit%covariance(ncoef,ncoef))
      fit%coefficients=coef; fit%conditional_variance=h; fit%residuals=resid
      fit%negative_log_likelihood=obj; fit%iterations=it; fit%status=status
      if(status==0) then
         fit%message='converged'
      else
         fit%message='optimizer reached iteration limit'
      end if
      call numerical_hessian(objective,theta,hess,status)
      if(status==0) call invert_matrix(hess,hinv,status)
      if(status==0) then
         fit%covariance=hinv
      else
         fit%covariance=0.0_dp
      end if

   contains
      function objective(th) result(value)
         real(dp), intent(in) :: th(:)
         real(dp) :: value
         real(dp) :: c(size(th)),hv(size(x))
         call theta_to_coefficients(th,c)
         call garch_variance(x,c,p,q,hv)
         value=0.5_dp*sum(log(hv(maxpq+1:))+x(maxpq+1:)**2/hv(maxpq+1:))
         if(.not.(value<huge(value))) value=huge(value)/100.0_dp
      end function objective
   end function garch_fit

   subroutine theta_to_coefficients(theta,coef)
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: coef(:)
      real(dp), allocatable :: e(:)
      real(dp) :: denom
      integer :: n
      n=size(theta)
      coef(1)=exp(max(-50.0_dp,min(50.0_dp,theta(1))))
      if(n>1) then
         allocate(e(n-1))
         e=exp(max(-50.0_dp,min(50.0_dp,theta(2:))))
         denom=1.0_dp+sum(e)
         coef(2:)=0.999_dp*e/denom
      end if
   end subroutine theta_to_coefficients

   subroutine coefficients_to_theta(coef,theta)
      real(dp), intent(in) :: coef(:)
      real(dp), intent(out) :: theta(:)
      real(dp) :: remaining
      theta(1)=log(max(coef(1),1.0e-12_dp))
      if(size(coef)>1) then
         remaining=max(1.0e-8_dp,0.999_dp-sum(coef(2:)))
         theta(2:)=log(max(coef(2:),1.0e-8_dp)/remaining)
      end if
   end subroutine coefficients_to_theta

end module tseries_models
