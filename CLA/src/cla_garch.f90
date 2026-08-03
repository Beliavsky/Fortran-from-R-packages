! SPDX-License-Identifier: GPL-3.0-or-later
module cla_garch
   use kind_mod, only: dp
   use cla_types, only: cla_garch_result_t, cla_success, cla_invalid_input, cla_garch_failure
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: cla_distribution_normal = 1
   integer, parameter, public :: cla_distribution_student = 2
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: cla_mu_sigma_garch

contains

   function cla_mu_sigma_garch(prices, arch_order, garch_order, distribution, &
      max_iterations, tolerance) result(out)
      !! Forecast asset means and covariance using separate constant-mean GARCH models.
      real(dp), intent(in) :: prices(:,:)
      integer, intent(in), optional :: arch_order, garch_order, distribution
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(cla_garch_result_t) :: out

      real(dp), allocatable :: returns(:,:), standardized(:,:), sigma(:), residual(:)
      real(dp), allocatable :: alpha(:), beta(:)
      integer :: nobs, nassets, p, q, dist, limit, asset, i, j, status
      real(dp) :: tol, corr, mean_i, mean_j, sd_i, sd_j, location, omega, shape
      real(dp) :: forecast_variance
      logical :: converged

      nobs = size(prices,1)
      nassets = size(prices,2)
      p = 1
      q = 1
      dist = cla_distribution_student
      limit = 600
      tol = 1.0e-6_dp
      if (present(arch_order)) p = arch_order
      if (present(garch_order)) q = garch_order
      if (present(distribution)) dist = distribution
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) tol = tolerance

      if (nobs < 21 .or. nassets < 2 .or. p < 0 .or. q < 0 .or. max(p,q) < 1 .or. &
          limit < 1 .or. tol <= 0.0_dp .or. any(prices <= 0.0_dp) .or. &
          .not. all(ieee_is_finite(prices)) .or. &
          (dist /= cla_distribution_normal .and. dist /= cla_distribution_student)) then
         out%info = cla_invalid_input
         return
      end if

      allocate(returns(nobs-1,nassets), standardized(nobs-1,nassets))
      do asset = 1, nassets
         returns(:,asset) = log(prices(2:nobs,asset)/prices(1:nobs-1,asset))
      end do
      allocate(out%mu(nassets), out%covariance(nassets,nassets), &
         out%forecast_sigma(nassets), out%conditional_sigma(nobs-1,nassets), &
         out%fitted_mean(nassets), out%omega(nassets), out%shape(nassets), &
         out%alpha(p,nassets), out%beta(q,nassets), out%converged(nassets))
      out%mu = 0.0_dp
      out%covariance = 0.0_dp
      out%forecast_sigma = 0.0_dp
      out%conditional_sigma = 0.0_dp
      out%fitted_mean = 0.0_dp
      out%omega = 0.0_dp
      out%shape = 0.0_dp
      out%alpha = 0.0_dp
      out%beta = 0.0_dp
      out%converged = .false.

      do asset = 1, nassets
         call fit_garch(-returns(:,asset),p,q,dist,limit,tol,location,omega,alpha,beta, &
            shape,sigma,residual,forecast_variance,converged,status)
         if (status /= 0) then
            out%info = cla_garch_failure
            return
         end if
         out%mu(asset) = -location
         out%forecast_sigma(asset) = sqrt(forecast_variance)
         out%conditional_sigma(:,asset) = sigma
         out%fitted_mean(asset) = location
         out%omega(asset) = omega
         out%shape(asset) = shape
         if (p > 0) out%alpha(:,asset) = alpha
         if (q > 0) out%beta(:,asset) = beta
         out%converged(asset) = converged
         standardized(:,asset) = returns(:,asset)/max(sigma,sqrt(tiny(1.0_dp)))
      end do

      do j = 1, nassets
         mean_j = sum(standardized(:,j))/real(nobs-1,dp)
         sd_j = sqrt(max(sum((standardized(:,j)-mean_j)**2)/real(nobs-2,dp),0.0_dp))
         do i = 1, j
            mean_i = sum(standardized(:,i))/real(nobs-1,dp)
            sd_i = sqrt(max(sum((standardized(:,i)-mean_i)**2)/real(nobs-2,dp),0.0_dp))
            if (sd_i > 0.0_dp .and. sd_j > 0.0_dp) then
               corr = sum((standardized(:,i)-mean_i)*(standardized(:,j)-mean_j))/ &
                  (real(nobs-2,dp)*sd_i*sd_j)
            else
               corr = merge(1.0_dp,0.0_dp,i==j)
            end if
            out%covariance(i,j) = out%forecast_sigma(i)*corr*out%forecast_sigma(j)
            out%covariance(j,i) = out%covariance(i,j)
         end do
      end do
      out%info = cla_success
   end function cla_mu_sigma_garch

   subroutine fit_garch(series,p,q,distribution,max_iterations,tolerance,location,omega, &
      alpha,beta,shape,sigma,residual,forecast_variance,converged,info)
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: p,q,distribution,max_iterations
      real(dp), intent(in) :: tolerance
      real(dp), intent(out) :: location,omega,shape,forecast_variance
      real(dp), allocatable, intent(out) :: alpha(:),beta(:),sigma(:),residual(:)
      logical, intent(out) :: converged
      integer, intent(out) :: info
      real(dp), allocatable :: start(:), optimum(:), variance_path(:)
      real(dp) :: sample_mean, sample_variance, total, remainder
      integer :: d, i

      sample_mean = sum(series)/real(size(series),dp)
      sample_variance = sum((series-sample_mean)**2)/real(max(1,size(series)-1),dp)
      sample_variance = max(sample_variance,1.0e-10_dp)
      d = 2+p+q+merge(1,0,distribution==cla_distribution_student)
      allocate(start(d))
      start = 0.0_dp
      start(1) = sample_mean
      start(2) = log(max(0.05_dp*sample_variance,1.0e-12_dp))
      total = 0.0_dp
      if (p > 0) total = total + 0.10_dp
      if (q > 0) total = total + 0.80_dp
      remainder = 0.995_dp-total
      if (p > 0) then
         do i=1,p
            start(2+i)=log((0.10_dp/real(p,dp))/remainder)
         end do
      end if
      if (q > 0) then
         do i=1,q
            start(2+p+i)=log((0.80_dp/real(q,dp))/remainder)
         end do
      end if
      if (distribution==cla_distribution_student) start(d)=log(8.0_dp-2.05_dp)

      call nelder_mead_garch(series,p,q,distribution,start,max_iterations,tolerance, &
         optimum,converged,info)
      if (info /= 0) return
      call unpack_parameters(optimum,p,q,distribution,location,omega,alpha,beta,shape)
      call garch_filter(series,location,omega,alpha,beta,shape,distribution, &
         variance_path,residual,forecast_variance,info)
      if (info /= 0) return
      allocate(sigma(size(series)))
      sigma = sqrt(variance_path)
   end subroutine fit_garch

   subroutine nelder_mead_garch(series,p,q,distribution,start,max_iterations,tolerance, &
      optimum,converged,info)
      real(dp), intent(in) :: series(:),start(:),tolerance
      integer, intent(in) :: p,q,distribution,max_iterations
      real(dp), allocatable, intent(out) :: optimum(:)
      logical, intent(out) :: converged
      integer, intent(out) :: info
      real(dp), allocatable :: simplex(:,:), values(:), centroid(:), reflected(:), &
         expanded(:), contracted(:), tmp(:)
      real(dp) :: fr,fe,fc,spread,step
      integer :: d,j,iteration,best,worst,second
      real(dp), parameter :: reflection=1.0_dp, expansion=2.0_dp, &
         contraction=0.5_dp, shrink=0.5_dp

      d=size(start)
      allocate(simplex(d,d+1),values(d+1),centroid(d),reflected(d),expanded(d), &
         contracted(d),tmp(d))
      simplex(:,1)=start
      do j=1,d
         simplex(:,j+1)=start
         step=0.15_dp*max(1.0_dp,abs(start(j)))
         simplex(j,j+1)=simplex(j,j+1)+step
      end do
      do j=1,d+1
         values(j)=garch_objective(series,p,q,distribution,simplex(:,j))
      end do
      converged=.false.
      do iteration=1,max_iterations
         call sort_simplex(simplex,values)
         spread=maxval(abs(simplex(:,2:d+1)-spread_matrix(simplex(:,1),d)))
         if(maxval(abs(values(2:d+1)-values(1)))<=tolerance .and. &
            spread<=sqrt(tolerance))then
            converged=.true.
            exit
         end if
         best=1; worst=d+1; second=d
         centroid=sum(simplex(:,1:d),dim=2)/real(d,dp)
         reflected=centroid+reflection*(centroid-simplex(:,worst))
         fr=garch_objective(series,p,q,distribution,reflected)
         if(fr<values(best))then
            expanded=centroid+expansion*(reflected-centroid)
            fe=garch_objective(series,p,q,distribution,expanded)
            if(fe<fr)then
               simplex(:,worst)=expanded;values(worst)=fe
            else
               simplex(:,worst)=reflected;values(worst)=fr
            end if
         else if(fr<values(second))then
            simplex(:,worst)=reflected;values(worst)=fr
         else
            if(fr<values(worst))then
               contracted=centroid+contraction*(reflected-centroid)
            else
               contracted=centroid+contraction*(simplex(:,worst)-centroid)
            end if
            fc=garch_objective(series,p,q,distribution,contracted)
            if(fc<min(fr,values(worst)))then
               simplex(:,worst)=contracted;values(worst)=fc
            else
               do j=2,d+1
                  simplex(:,j)=simplex(:,best)+shrink*(simplex(:,j)-simplex(:,best))
                  values(j)=garch_objective(series,p,q,distribution,simplex(:,j))
               end do
            end if
         end if
      end do
      call sort_simplex(simplex,values)
      allocate(optimum(d));optimum=simplex(:,1)
      info=merge(0,1,ieee_is_finite(values(1)) .and. values(1)<0.5_dp*huge(1.0_dp))
   contains
      pure function spread_matrix(column,ncol) result(matrix)
         real(dp),intent(in)::column(:)
         integer,intent(in)::ncol
         real(dp)::matrix(size(column),ncol)
         integer::k
         do k=1,ncol
            matrix(:,k)=column
         end do
      end function spread_matrix
      subroutine sort_simplex(x,f)
         real(dp),intent(inout)::x(:,:),f(:)
         integer::a,b
         real(dp)::fv
         do a=2,size(f)
            tmp=x(:,a);fv=f(a);b=a-1
            do while(b>=1)
               if(f(b)<=fv)exit
               x(:,b+1)=x(:,b);f(b+1)=f(b);b=b-1
            end do
            x(:,b+1)=tmp;f(b+1)=fv
         end do
      end subroutine sort_simplex
   end subroutine nelder_mead_garch

   real(dp) function garch_objective(series,p,q,distribution,raw) result(value)
      real(dp),intent(in)::series(:),raw(:)
      integer,intent(in)::p,q,distribution
      real(dp)::location,omega,shape,forecast
      real(dp),allocatable::alpha(:),beta(:),variance_path(:),residual(:)
      integer::status,t
      call unpack_parameters(raw,p,q,distribution,location,omega,alpha,beta,shape)
      call garch_filter(series,location,omega,alpha,beta,shape,distribution, &
         variance_path,residual,forecast,status)
      if(status/=0)then
         value=huge(1.0_dp)
         return
      end if
      value=0.0_dp
      do t=1,size(series)
         value=value-log_density(residual(t),variance_path(t),shape,distribution)
      end do
      if(.not.ieee_is_finite(value))value=huge(1.0_dp)
   end function garch_objective

   subroutine unpack_parameters(raw,p,q,distribution,location,omega,alpha,beta,shape)
      real(dp),intent(in)::raw(:)
      integer,intent(in)::p,q,distribution
      real(dp),intent(out)::location,omega,shape
      real(dp),allocatable,intent(out)::alpha(:),beta(:)
      real(dp),allocatable::e(:)
      real(dp)::denominator
      integer::k
      location=raw(1)
      omega=exp(max(-40.0_dp,min(20.0_dp,raw(2))))
      allocate(alpha(p),beta(q),e(p+q))
      do k=1,p+q
         e(k)=exp(max(-30.0_dp,min(30.0_dp,raw(2+k))))
      end do
      denominator=1.0_dp+sum(e)
      if(p>0)alpha=0.995_dp*e(1:p)/denominator
      if(q>0)beta=0.995_dp*e(p+1:p+q)/denominator
      shape=0.0_dp
      if(distribution==cla_distribution_student) &
         shape=2.05_dp+exp(max(-20.0_dp,min(5.0_dp,raw(size(raw)))))
   end subroutine unpack_parameters

   subroutine garch_filter(series,location,omega,alpha,beta,shape,distribution, &
      variance_path,residual,forecast_variance,info)
      real(dp),intent(in)::series(:),location,omega,alpha(:),beta(:),shape
      integer,intent(in)::distribution
      real(dp),allocatable,intent(out)::variance_path(:),residual(:)
      real(dp),intent(out)::forecast_variance
      integer,intent(out)::info
      real(dp)::initial_variance,value
      integer::t,lag,n,p,q
      n=size(series);p=size(alpha);q=size(beta)
      allocate(variance_path(n),residual(n))
      residual=series-location
      initial_variance=max(sum(residual**2)/real(max(1,n),dp),1.0e-12_dp)
      do t=1,n
         value=omega
         do lag=1,p
            if(t-lag>=1)then
               value=value+alpha(lag)*residual(t-lag)**2
            else
               value=value+alpha(lag)*initial_variance
            end if
         end do
         do lag=1,q
            if(t-lag>=1)then
               value=value+beta(lag)*variance_path(t-lag)
            else
               value=value+beta(lag)*initial_variance
            end if
         end do
         variance_path(t)=max(value,1.0e-14_dp)
      end do
      forecast_variance=omega
      do lag=1,p
         if(n-lag+1>=1)then
            forecast_variance=forecast_variance+alpha(lag)*residual(n-lag+1)**2
         else
            forecast_variance=forecast_variance+alpha(lag)*initial_variance
         end if
      end do
      do lag=1,q
         if(n-lag+1>=1)then
            forecast_variance=forecast_variance+beta(lag)*variance_path(n-lag+1)
         else
            forecast_variance=forecast_variance+beta(lag)*initial_variance
         end if
      end do
      forecast_variance=max(forecast_variance,1.0e-14_dp)
      info=0
      if(.not.all(ieee_is_finite(variance_path)) .or. &
         .not.ieee_is_finite(forecast_variance) .or. omega<=0.0_dp .or. &
         (distribution==cla_distribution_student .and. shape<=2.0_dp))info=1
   end subroutine garch_filter

   real(dp) function log_density(residual,variance_value,shape,distribution) result(value)
      real(dp),intent(in)::residual,variance_value,shape
      integer,intent(in)::distribution
      real(dp)::z
      z=residual/sqrt(variance_value)
      select case(distribution)
      case(cla_distribution_normal)
         value=-0.5_dp*(log(2.0_dp*pi)+log(variance_value)+z*z)
      case(cla_distribution_student)
         value=log_gamma(0.5_dp*(shape+1.0_dp))-log_gamma(0.5_dp*shape) &
            -0.5_dp*log((shape-2.0_dp)*pi)-0.5_dp*log(variance_value) &
            -0.5_dp*(shape+1.0_dp)*log(1.0_dp+z*z/(shape-2.0_dp))
      case default
         value=-huge(1.0_dp)
      end select
   end function log_density

end module cla_garch
