! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_gogarch
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : normalize_covariance
   use rmgarch_rng, only : random_normal
   use rmgarch_types, only : gogarch_fit_result
   use rmgarch_univariate, only : fit_marginal_garch11, filter_garch11
   use rmgarch_ica, only : fastica, radical
   implicit none
   private

   public :: gogarch_covariance, gogarch_correlation, gogarch_sigma
   public :: gogarch_coskewness, gogarch_cokurtosis
   public :: portfolio_factor_moments, portfolio_covariance_beta
   public :: portfolio_coskew_beta, portfolio_cokurt_beta
   public :: simulate_gogarch, gogarch_moments_at
   public :: fit_gogarch11, filter_gogarch11, forecast_gogarch11, simulate_fitted_gogarch11

contains

   subroutine gogarch_covariance(component_sigma, mixing, covariance)
      real(dp), intent(in) :: component_sigma(:,:), mixing(:,:)
      real(dp), intent(out) :: covariance(size(mixing,1),size(mixing,1),size(component_sigma,1))
      real(dp) :: scaled(size(mixing,1),size(mixing,2))
      integer :: t, j
      do t = 1, size(component_sigma,1)
         scaled = mixing
         do j = 1, size(mixing,2)
            scaled(:,j) = scaled(:,j)*component_sigma(t,j)**2
         end do
         covariance(:,:,t) = matmul(scaled,transpose(mixing))
      end do
   end subroutine gogarch_covariance

   subroutine gogarch_correlation(component_sigma, mixing, correlation)
      real(dp), intent(in) :: component_sigma(:,:), mixing(:,:)
      real(dp), intent(out) :: correlation(size(mixing,1),size(mixing,1),size(component_sigma,1))
      real(dp), allocatable :: covariance(:,:,:)
      integer :: t
      allocate(covariance(size(mixing,1),size(mixing,1),size(component_sigma,1)))
      call gogarch_covariance(component_sigma,mixing,covariance)
      do t = 1, size(component_sigma,1)
         correlation(:,:,t) = normalize_covariance(covariance(:,:,t))
      end do
   end subroutine gogarch_correlation

   subroutine gogarch_sigma(component_sigma, mixing, sigma)
      real(dp), intent(in) :: component_sigma(:,:), mixing(:,:)
      real(dp), intent(out) :: sigma(size(component_sigma,1),size(mixing,1))
      real(dp), allocatable :: covariance(:,:,:)
      integer :: t, i
      allocate(covariance(size(mixing,1),size(mixing,1),size(component_sigma,1)))
      call gogarch_covariance(component_sigma,mixing,covariance)
      do t = 1, size(component_sigma,1)
         do i = 1, size(mixing,1)
            sigma(t,i) = sqrt(max(covariance(i,i,t),0.0_dp))
         end do
      end do
   end subroutine gogarch_sigma

   function gogarch_coskewness(mixing, factor_third_moment) result(coskew)
      real(dp), intent(in) :: mixing(:,:), factor_third_moment(:)
      real(dp) :: coskew(size(mixing,1),size(mixing,1),size(mixing,1))
      integer :: i, j, k, f, m
      m = size(mixing,1); coskew = 0.0_dp
      do f = 1, min(size(mixing,2),size(factor_third_moment))
         do k = 1, m
            do j = 1, m
               do i = 1, m
                  coskew(i,j,k) = coskew(i,j,k)+factor_third_moment(f)* &
                     mixing(i,f)*mixing(j,f)*mixing(k,f)
               end do
            end do
         end do
      end do
   end function gogarch_coskewness

   function gogarch_cokurtosis(mixing, factor_variance, factor_excess_kurtosis) result(cokurt)
      real(dp), intent(in) :: mixing(:,:), factor_variance(:), factor_excess_kurtosis(:)
      real(dp) :: cokurt(size(mixing,1),size(mixing,1),size(mixing,1),size(mixing,1))
      real(dp) :: cov(size(mixing,1),size(mixing,1)), cumulant
      integer :: i, j, k, l, f, m
      m = size(mixing,1); cov = 0.0_dp; cokurt = 0.0_dp
      do f = 1, min(size(mixing,2),size(factor_variance))
         do j = 1, m
            do i = 1, m
               cov(i,j) = cov(i,j)+mixing(i,f)*mixing(j,f)*factor_variance(f)
            end do
         end do
      end do
      do l = 1, m
         do k = 1, m
            do j = 1, m
               do i = 1, m
                  cokurt(i,j,k,l) = cov(i,j)*cov(k,l)+cov(i,k)*cov(j,l)+cov(i,l)*cov(j,k)
               end do
            end do
         end do
      end do
      do f = 1, min(size(mixing,2),min(size(factor_variance),size(factor_excess_kurtosis)))
         cumulant = factor_excess_kurtosis(f)*factor_variance(f)**2
         do l = 1, m
            do k = 1, m
               do j = 1, m
                  do i = 1, m
                     cokurt(i,j,k,l) = cokurt(i,j,k,l)+cumulant*mixing(i,f)*mixing(j,f)*mixing(k,f)*mixing(l,f)
                  end do
               end do
            end do
         end do
      end do
   end function gogarch_cokurtosis

   subroutine portfolio_factor_moments(weights, mixing, factor_mean, factor_variance, factor_third, &
      factor_excess_kurtosis, mean, variance, skewness, kurtosis)
      real(dp), intent(in) :: weights(:), mixing(:,:), factor_mean(:), factor_variance(:)
      real(dp), intent(in) :: factor_third(:), factor_excess_kurtosis(:)
      real(dp), intent(out) :: mean, variance, skewness, kurtosis
      real(dp) :: exposure(size(mixing,2)), third, fourth_cumulant
      exposure = matmul(transpose(mixing),weights)
      mean = sum(exposure*factor_mean)
      variance = sum(exposure**2*factor_variance)
      third = sum(exposure**3*factor_third)
      fourth_cumulant = sum(exposure**4*factor_excess_kurtosis*factor_variance**2)
      if (variance > 0.0_dp) then
         skewness = third/variance**1.5_dp
         kurtosis = 3.0_dp+fourth_cumulant/variance**2
      else
         skewness = 0.0_dp; kurtosis = 0.0_dp
      end if
   end subroutine portfolio_factor_moments

   function portfolio_covariance_beta(weights, covariance) result(beta)
      real(dp), intent(in) :: weights(:), covariance(:,:)
      real(dp) :: beta(size(weights)), portfolio_variance
      portfolio_variance = dot_product(weights,matmul(covariance,weights))
      if (portfolio_variance > 0.0_dp) then
         beta = matmul(covariance,weights)/portfolio_variance
      else
         beta = 0.0_dp
      end if
   end function portfolio_covariance_beta

   function portfolio_coskew_beta(weights, coskew) result(beta)
      !! Co-skewness beta for each asset relative to the weighted portfolio.
      real(dp), intent(in) :: weights(:)
      real(dp), intent(in) :: coskew(:,:,:)
      real(dp) :: beta(size(weights)), numerator(size(weights)), denominator
      integer :: i, j, k, m

      m = size(weights)
      numerator = 0.0_dp
      if (size(coskew,1) /= m .or. size(coskew,2) /= m .or. size(coskew,3) /= m) then
         beta = 0.0_dp
         return
      end if
      do i = 1, m
         do k = 1, m
            do j = 1, m
               numerator(i) = numerator(i)+coskew(i,j,k)*weights(j)*weights(k)
            end do
         end do
      end do
      denominator = dot_product(weights,numerator)
      if (abs(denominator) > tiny(1.0_dp)) then
         beta = numerator/denominator
      else
         beta = 0.0_dp
      end if
   end function portfolio_coskew_beta

   function portfolio_cokurt_beta(weights, cokurt) result(beta)
      !! Co-kurtosis beta for each asset relative to the weighted portfolio.
      real(dp), intent(in) :: weights(:)
      real(dp), intent(in) :: cokurt(:,:,:,:)
      real(dp) :: beta(size(weights)), numerator(size(weights)), denominator
      integer :: i, j, k, l, m

      m = size(weights)
      numerator = 0.0_dp
      if (size(cokurt,1) /= m .or. size(cokurt,2) /= m .or. &
          size(cokurt,3) /= m .or. size(cokurt,4) /= m) then
         beta = 0.0_dp
         return
      end if
      do i = 1, m
         do l = 1, m
            do k = 1, m
               do j = 1, m
                  numerator(i) = numerator(i)+cokurt(i,j,k,l)*weights(j)*weights(k)*weights(l)
               end do
            end do
         end do
      end do
      denominator = dot_product(weights,numerator)
      if (abs(denominator) > tiny(1.0_dp)) then
         beta = numerator/denominator
      else
         beta = 0.0_dp
      end if
   end function portfolio_cokurt_beta

   subroutine simulate_gogarch(component_shocks, component_mean, component_sigma, mixing, asset_returns, factors)
      !! Mix independent standardized factor shocks into asset returns.
      real(dp), intent(in) :: component_shocks(:,:), component_mean(:,:), component_sigma(:,:), mixing(:,:)
      real(dp), intent(out) :: asset_returns(size(component_shocks,1),size(mixing,1))
      real(dp), intent(out), optional :: factors(size(component_shocks,1),size(component_shocks,2))
      real(dp) :: factor_values(size(component_shocks,1),size(component_shocks,2))

      if (size(component_mean,1) /= size(component_shocks,1) .or. &
          size(component_mean,2) /= size(component_shocks,2) .or. &
          size(component_sigma,1) /= size(component_shocks,1) .or. &
          size(component_sigma,2) /= size(component_shocks,2) .or. &
          size(mixing,2) /= size(component_shocks,2)) then
         asset_returns = 0.0_dp
         if (present(factors)) factors = 0.0_dp
         return
      end if
      factor_values = component_mean+component_sigma*component_shocks
      asset_returns = matmul(factor_values,transpose(mixing))
      if (present(factors)) factors = factor_values
   end subroutine simulate_gogarch

   subroutine gogarch_moments_at(mixing, factor_mean, factor_variance, factor_third, &
      factor_excess_kurtosis, asset_mean, covariance, coskew, cokurt)
      !! Construct one-period GO-GARCH first through fourth co-moments.
      real(dp), intent(in) :: mixing(:,:), factor_mean(:), factor_variance(:)
      real(dp), intent(in) :: factor_third(:), factor_excess_kurtosis(:)
      real(dp), intent(out) :: asset_mean(size(mixing,1))
      real(dp), intent(out) :: covariance(size(mixing,1),size(mixing,1))
      real(dp), intent(out) :: coskew(size(mixing,1),size(mixing,1),size(mixing,1))
      real(dp), intent(out) :: cokurt(size(mixing,1),size(mixing,1),size(mixing,1),size(mixing,1))
      real(dp) :: component_sigma(1,size(mixing,2)), cov3(size(mixing,1),size(mixing,1),1)

      asset_mean = matmul(mixing,factor_mean)
      component_sigma(1,:) = sqrt(max(factor_variance,0.0_dp))
      call gogarch_covariance(component_sigma,mixing,cov3)
      covariance = cov3(:,:,1)
      coskew = gogarch_coskewness(mixing,factor_third)
      cokurt = gogarch_cokurtosis(mixing,factor_variance,factor_excess_kurtosis)
   end subroutine gogarch_moments_at

   function fit_gogarch11(data, ica_method, max_iterations) result(fit)
      !! Fit a square GO-GARCH approximation using ICA followed by independent
      !! Gaussian GARCH(1,1) component models.
      real(dp), intent(in) :: data(:,:)
      character(len=*), intent(in), optional :: ica_method
      integer, intent(in), optional :: max_iterations
      type(gogarch_fit_result) :: fit
      character(len=16) :: method
      integer :: n, m, j, maxit

      n = size(data,1)
      m = size(data,2)
      maxit = 500
      if (present(max_iterations)) maxit = max(1,max_iterations)
      method = 'fastica'
      if (present(ica_method)) method = lowercase(adjustl(ica_method))
      if (n <= 2 .or. m < 1) then
         fit%status = 2
         return
      end if
      select case (trim(method))
      case ('radical')
         fit%ica = radical(data,max_sweeps=8,angle_points=41)
      case default
         fit%ica = fastica(data,max_iterations=maxit)
      end select
      if (fit%ica%status > 1) then
         fit%status = fit%ica%status
         return
      end if
      allocate(fit%components(m),fit%component_sigma(n,m),fit%standardized_components(n,m))
      call fit_marginal_garch11(fit%ica%sources,fit%components,fit%standardized_components, &
         fit%component_sigma,maxit)
      allocate(fit%covariance(m,m,n),fit%correlation(m,m,n))
      call gogarch_covariance(fit%component_sigma,fit%ica%mixing,fit%covariance)
      call gogarch_correlation(fit%component_sigma,fit%ica%mixing,fit%correlation)
      fit%status = fit%ica%status
      do j = 1, m
         fit%status = max(fit%status,fit%components(j)%status)
      end do
   end function fit_gogarch11

   subroutine filter_gogarch11(data, fit, component_returns, component_sigma, &
      standardized_components, covariance, correlation, valid)
      !! Apply fitted ICA and component GARCH(1,1) parameters to another data set.
      real(dp), intent(in) :: data(:,:)
      type(gogarch_fit_result), intent(in) :: fit
      real(dp), intent(out) :: component_returns(size(data,1),size(data,2))
      real(dp), intent(out) :: component_sigma(size(data,1),size(data,2))
      real(dp), intent(out) :: standardized_components(size(data,1),size(data,2))
      real(dp), intent(out) :: covariance(size(data,2),size(data,2),size(data,1))
      real(dp), intent(out) :: correlation(size(data,2),size(data,2),size(data,1))
      logical, intent(out), optional :: valid
      real(dp), allocatable :: residuals(:)
      real(dp) :: log_likelihood
      integer :: n, m, j
      logical :: ok, margin_ok

      n = size(data,1)
      m = size(data,2)
      component_returns = 0.0_dp
      component_sigma = 0.0_dp
      standardized_components = 0.0_dp
      covariance = 0.0_dp
      correlation = 0.0_dp
      ok = n >= 2 .and. allocated(fit%components) .and. size(fit%components) == m .and. &
         allocated(fit%ica%center) .and. size(fit%ica%center) == m .and. &
         allocated(fit%ica%unmixing) .and. all(shape(fit%ica%unmixing) == [m,m]) .and. &
         allocated(fit%ica%mixing) .and. all(shape(fit%ica%mixing) == [m,m])
      if (.not. ok) then
         if (present(valid)) valid = .false.
         return
      end if

      component_returns = matmul(data-spread(fit%ica%center,1,n),transpose(fit%ica%unmixing))
      allocate(residuals(n))
      do j = 1, m
         call filter_garch11(component_returns(:,j),fit%components(j)%mean, &
            fit%components(j)%omega,fit%components(j)%alpha,fit%components(j)%beta, &
            residuals,component_sigma(:,j),standardized_components(:,j), &
            log_likelihood,margin_ok)
         ok = ok .and. margin_ok
      end do
      call gogarch_covariance(component_sigma,fit%ica%mixing,covariance)
      call gogarch_correlation(component_sigma,fit%ica%mixing,correlation)
      if (present(valid)) valid = ok
   end subroutine filter_gogarch11

   subroutine forecast_gogarch11(fit, horizons, component_sigma, covariance, correlation)
      type(gogarch_fit_result), intent(in) :: fit
      integer, intent(in) :: horizons
      real(dp), intent(out) :: component_sigma(horizons,size(fit%components))
      real(dp), intent(out) :: covariance(size(fit%ica%mixing,1),size(fit%ica%mixing,1),horizons)
      real(dp), intent(out) :: correlation(size(fit%ica%mixing,1),size(fit%ica%mixing,1),horizons)
      real(dp) :: variance
      integer :: h, j

      component_sigma = 0.0_dp
      covariance = 0.0_dp
      correlation = 0.0_dp
      if (horizons <= 0 .or. .not. allocated(fit%components)) return
      do j = 1, size(fit%components)
         if (.not. allocated(fit%components(j)%residuals) .or. &
             .not. allocated(fit%components(j)%sigma)) cycle
         variance = fit%components(j)%omega+fit%components(j)%alpha* &
            fit%components(j)%residuals(size(fit%components(j)%residuals))**2+ &
            fit%components(j)%beta*fit%components(j)%sigma(size(fit%components(j)%sigma))**2
         component_sigma(1,j) = sqrt(max(variance,1.0e-12_dp))
         do h = 2, horizons
            variance = fit%components(j)%omega+(fit%components(j)%alpha+ &
               fit%components(j)%beta)*variance
            component_sigma(h,j) = sqrt(max(variance,1.0e-12_dp))
         end do
      end do
      call gogarch_covariance(component_sigma,fit%ica%mixing,covariance)
      call gogarch_correlation(component_sigma,fit%ica%mixing,correlation)
   end subroutine forecast_gogarch11

   subroutine simulate_fitted_gogarch11(nobs, fit, asset_returns, component_returns, burn)
      integer, intent(in) :: nobs
      type(gogarch_fit_result), intent(in) :: fit
      real(dp), intent(out) :: asset_returns(nobs,size(fit%ica%mixing,1))
      real(dp), intent(out), optional :: component_returns(nobs,size(fit%components))
      integer, intent(in), optional :: burn
      real(dp), allocatable :: factors(:,:), h(:), residual(:)
      real(dp) :: innovation
      integer :: nburn, ntotal, t, j, m

      asset_returns = 0.0_dp
      if (nobs < 1 .or. .not. allocated(fit%components) .or. &
          .not. allocated(fit%ica%mixing) .or. .not. allocated(fit%ica%center)) then
         if (present(component_returns)) component_returns = 0.0_dp
         return
      end if
      nburn = 300
      if (present(burn)) nburn = max(0,burn)
      ntotal = nobs+nburn
      m = size(fit%components)
      allocate(factors(ntotal,m),h(m),residual(m))
      do j = 1, m
         h(j) = fit%components(j)%omega/max(1.0_dp-fit%components(j)%alpha- &
            fit%components(j)%beta,1.0e-6_dp)
         residual(j) = 0.0_dp
      end do
      do t = 1, ntotal
         do j = 1, m
            if (t > 1) h(j) = fit%components(j)%omega+fit%components(j)%alpha* &
               residual(j)**2+fit%components(j)%beta*h(j)
            innovation = random_normal()
            residual(j) = sqrt(max(h(j),1.0e-12_dp))*innovation
            factors(t,j) = fit%components(j)%mean+residual(j)
         end do
      end do
      asset_returns = matmul(factors(nburn+1:ntotal,:),transpose(fit%ica%mixing))+ &
         spread(fit%ica%center,1,nobs)
      if (present(component_returns)) component_returns = factors(nburn+1:ntotal,:)
   end subroutine simulate_fitted_gogarch11

   pure function lowercase(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, code
      out = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code+32)
      end do
   end function lowercase

end module rmgarch_gogarch
