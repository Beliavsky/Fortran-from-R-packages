! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_rolling
   use rmgarch_kinds, only : dp
   use rmgarch_types, only : dcc_fit_result, multivariate_garch_fit, &
      gogarch_fit_result, rolling_dcc_result, rolling_gogarch_result, dist_gaussian
   use rmgarch_dcc, only : fit_dcc, dcc_forecast_history
   use rmgarch_model, only : fit_two_step_dcc_general
   use rmgarch_gogarch, only : fit_gogarch11, forecast_gogarch11
   implicit none
   private

   public :: roll_dcc, roll_two_step_dcc, roll_gogarch11

contains

   function roll_dcc(z, window, refit_every, horizon, p, q, g, distribution, &
      shape, estimate_shape, max_iterations) result(rolling)
      !! Rolling DCC estimation for already-standardized score series.
      real(dp), intent(in) :: z(:,:)
      integer, intent(in) :: window, refit_every, horizon
      integer, intent(in), optional :: p, q, g, distribution, max_iterations
      real(dp), intent(in), optional :: shape
      logical, intent(in), optional :: estimate_shape
      type(rolling_dcc_result) :: rolling
      type(dcc_fit_result) :: fit
      real(dp), allocatable :: qforecast(:,:,:), rforecast(:,:,:)
      integer :: pp, qq, gg, dist, maxit, nfit, index, origin, m
      real(dp) :: nu
      logical :: fit_shape

      call rolling_options(p,q,g,distribution,shape,estimate_shape,max_iterations, &
         pp,qq,gg,dist,nu,fit_shape,maxit)
      rolling%window = window
      rolling%refit_every = refit_every
      rolling%horizon = horizon
      if (.not. valid_rolling_request(size(z,1),size(z,2),window,refit_every,horizon,pp,qq,gg)) then
         call allocate_empty_dcc_rolling(rolling,size(z,2),max(horizon,0),pp,qq,gg)
         return
      end if

      nfit = 1+(size(z,1)-window)/refit_every
      m = size(z,2)
      allocate(rolling%origin(nfit),rolling%fit_status(nfit))
      allocate(rolling%alpha(pp,nfit),rolling%beta(qq,nfit),rolling%gamma(gg,nfit))
      allocate(rolling%correlation(m,m,horizon,nfit))
      rolling%alpha = 0.0_dp
      rolling%beta = 0.0_dp
      rolling%gamma = 0.0_dp
      rolling%correlation = 0.0_dp

      do index = 1, nfit
         origin = window+(index-1)*refit_every
         rolling%origin(index) = origin
         fit = fit_dcc(z(origin-window+1:origin,:),p=pp,q=qq,g=gg,distribution=dist, &
            shape=nu,estimate_shape=fit_shape,max_iterations=maxit)
         rolling%fit_status(index) = fit%status
         if (allocated(fit%spec%alpha) .and. pp > 0) rolling%alpha(:,index) = fit%spec%alpha
         if (allocated(fit%spec%beta) .and. qq > 0) rolling%beta(:,index) = fit%spec%beta
         if (allocated(fit%spec%gamma) .and. gg > 0) rolling%gamma(:,index) = fit%spec%gamma
         if (allocated(fit%q) .and. allocated(fit%qbar) .and. allocated(fit%nbar)) then
            allocate(qforecast(m,m,horizon),rforecast(m,m,horizon))
            call dcc_forecast_history(fit%spec,fit%qbar,fit%nbar,fit%q, &
               z(origin-window+1:origin,:),horizon,qforecast,rforecast)
            rolling%correlation(:,:,:,index) = rforecast
            deallocate(qforecast,rforecast)
         end if
      end do
   end function roll_dcc

   function roll_two_step_dcc(data, window, refit_every, horizon, p, q, g, distribution, &
      shape, estimate_shape, max_iterations) result(rolling)
      !! Rolling raw-return workflow: Gaussian marginal GARCH(1,1), DCC fit,
      !! and multi-step conditional volatility/covariance/correlation forecasts.
      real(dp), intent(in) :: data(:,:)
      integer, intent(in) :: window, refit_every, horizon
      integer, intent(in), optional :: p, q, g, distribution, max_iterations
      real(dp), intent(in), optional :: shape
      logical, intent(in), optional :: estimate_shape
      type(rolling_dcc_result) :: rolling
      type(multivariate_garch_fit), allocatable :: fit
      real(dp), allocatable :: qforecast(:,:,:), rforecast(:,:,:)
      real(dp) :: variance
      integer :: pp, qq, gg, dist, maxit, nfit, index, origin, m, h, j, i
      real(dp) :: nu
      logical :: fit_shape

      call rolling_options(p,q,g,distribution,shape,estimate_shape,max_iterations, &
         pp,qq,gg,dist,nu,fit_shape,maxit)
      rolling%window = window
      rolling%refit_every = refit_every
      rolling%horizon = horizon
      if (.not. valid_rolling_request(size(data,1),size(data,2),window,refit_every,horizon,pp,qq,gg)) then
         call allocate_empty_dcc_rolling(rolling,size(data,2),max(horizon,0),pp,qq,gg)
         allocate(rolling%sigma(max(horizon,0),size(data,2),0))
         allocate(rolling%covariance(size(data,2),size(data,2),max(horizon,0),0))
         return
      end if

      nfit = 1+(size(data,1)-window)/refit_every
      m = size(data,2)
      allocate(fit)
      allocate(rolling%origin(nfit),rolling%fit_status(nfit))
      allocate(rolling%alpha(pp,nfit),rolling%beta(qq,nfit),rolling%gamma(gg,nfit))
      allocate(rolling%sigma(horizon,m,nfit),rolling%covariance(m,m,horizon,nfit))
      allocate(rolling%correlation(m,m,horizon,nfit))
      rolling%alpha = 0.0_dp
      rolling%beta = 0.0_dp
      rolling%gamma = 0.0_dp
      rolling%sigma = 0.0_dp
      rolling%covariance = 0.0_dp
      rolling%correlation = 0.0_dp

      do index = 1, nfit
         origin = window+(index-1)*refit_every
         rolling%origin(index) = origin
         fit = fit_two_step_dcc_general(data(origin-window+1:origin,:),pp,qq,gg, &
            distribution=dist,shape=nu,estimate_shape=fit_shape,max_iterations=maxit)
         rolling%fit_status(index) = fit%status
         if (allocated(fit%dcc%spec%alpha) .and. pp > 0) rolling%alpha(:,index) = fit%dcc%spec%alpha
         if (allocated(fit%dcc%spec%beta) .and. qq > 0) rolling%beta(:,index) = fit%dcc%spec%beta
         if (allocated(fit%dcc%spec%gamma) .and. gg > 0) rolling%gamma(:,index) = fit%dcc%spec%gamma

         do j = 1, m
            if (.not. allocated(fit%margins(j)%residuals) .or. &
                .not. allocated(fit%margins(j)%sigma)) cycle
            variance = fit%margins(j)%omega+fit%margins(j)%alpha* &
               fit%margins(j)%residuals(window)**2+fit%margins(j)%beta* &
               fit%margins(j)%sigma(window)**2
            rolling%sigma(1,j,index) = sqrt(max(variance,1.0e-12_dp))
            do h = 2, horizon
               variance = fit%margins(j)%omega+(fit%margins(j)%alpha+ &
                  fit%margins(j)%beta)*variance
               rolling%sigma(h,j,index) = sqrt(max(variance,1.0e-12_dp))
            end do
         end do

         if (allocated(fit%dcc%q) .and. allocated(fit%dcc%qbar) .and. allocated(fit%dcc%nbar)) then
            allocate(qforecast(m,m,horizon),rforecast(m,m,horizon))
            call dcc_forecast_history(fit%dcc%spec,fit%dcc%qbar,fit%dcc%nbar,fit%dcc%q, &
               fit%standardized,horizon,qforecast,rforecast)
            rolling%correlation(:,:,:,index) = rforecast
            do h = 1, horizon
               do j = 1, m
                  do i = 1, m
                     rolling%covariance(i,j,h,index) = rolling%sigma(h,i,index)* &
                        rforecast(i,j,h)*rolling%sigma(h,j,index)
                  end do
               end do
            end do
            deallocate(qforecast,rforecast)
         end if
      end do
   end function roll_two_step_dcc

   function roll_gogarch11(data, window, refit_every, horizon, ica_method, &
      max_iterations) result(rolling)
      !! Rolling square GO-GARCH approximation with ICA and Gaussian
      !! component GARCH(1,1) models.
      real(dp), intent(in) :: data(:,:)
      integer, intent(in) :: window, refit_every, horizon
      character(len=*), intent(in), optional :: ica_method
      integer, intent(in), optional :: max_iterations
      type(rolling_gogarch_result) :: rolling
      type(gogarch_fit_result) :: fit
      real(dp), allocatable :: component_sigma(:,:), covariance(:,:,:), correlation(:,:,:)
      character(len=16) :: method
      integer :: nfit, index, origin, m, maxit

      method = 'fastica'
      if (present(ica_method)) method = ica_method
      maxit = 500
      if (present(max_iterations)) maxit = max(1,max_iterations)
      rolling%window = window
      rolling%refit_every = refit_every
      rolling%horizon = horizon
      m = size(data,2)
      if (window < 3 .or. window > size(data,1) .or. refit_every < 1 .or. &
          horizon < 1 .or. m < 1) then
         allocate(rolling%origin(0),rolling%fit_status(0))
         allocate(rolling%covariance(m,m,max(horizon,0),0))
         allocate(rolling%correlation(m,m,max(horizon,0),0))
         return
      end if

      nfit = 1+(size(data,1)-window)/refit_every
      allocate(rolling%origin(nfit),rolling%fit_status(nfit))
      allocate(rolling%covariance(m,m,horizon,nfit),rolling%correlation(m,m,horizon,nfit))
      rolling%covariance = 0.0_dp
      rolling%correlation = 0.0_dp
      allocate(component_sigma(horizon,m),covariance(m,m,horizon),correlation(m,m,horizon))
      do index = 1, nfit
         origin = window+(index-1)*refit_every
         rolling%origin(index) = origin
         fit = fit_gogarch11(data(origin-window+1:origin,:),ica_method=trim(method), &
            max_iterations=maxit)
         rolling%fit_status(index) = fit%status
         if (allocated(fit%components) .and. allocated(fit%ica%mixing)) then
            call forecast_gogarch11(fit,horizon,component_sigma,covariance,correlation)
            rolling%covariance(:,:,:,index) = covariance
            rolling%correlation(:,:,:,index) = correlation
         end if
      end do
   end function roll_gogarch11

   subroutine rolling_options(p,q,g,distribution,shape,estimate_shape,max_iterations, &
      pp,qq,gg,dist,nu,fit_shape,maxit)
      integer, intent(in), optional :: p, q, g, distribution, max_iterations
      real(dp), intent(in), optional :: shape
      logical, intent(in), optional :: estimate_shape
      integer, intent(out) :: pp, qq, gg, dist, maxit
      real(dp), intent(out) :: nu
      logical, intent(out) :: fit_shape

      pp = 1
      qq = 1
      gg = 0
      dist = dist_gaussian
      nu = 8.0_dp
      fit_shape = .false.
      maxit = 500
      if (present(p)) pp = max(0,p)
      if (present(q)) qq = max(0,q)
      if (present(g)) gg = max(0,g)
      if (present(distribution)) dist = distribution
      if (present(shape)) nu = shape
      if (present(estimate_shape)) fit_shape = estimate_shape
      if (present(max_iterations)) maxit = max(1,max_iterations)
   end subroutine rolling_options

   pure logical function valid_rolling_request(n,m,window,refit_every,horizon,p,q,g)
      integer, intent(in) :: n, m, window, refit_every, horizon, p, q, g
      valid_rolling_request = window >= 3 .and. window <= n .and. refit_every >= 1 .and. &
         horizon >= 1 .and. p+q+g >= 1 .and. m >= 1
   end function valid_rolling_request

   subroutine allocate_empty_dcc_rolling(rolling,m,horizon,p,q,g)
      type(rolling_dcc_result), intent(inout) :: rolling
      integer, intent(in) :: m, horizon, p, q, g
      allocate(rolling%origin(0),rolling%fit_status(0))
      allocate(rolling%alpha(max(p,0),0),rolling%beta(max(q,0),0),rolling%gamma(max(g,0),0))
      allocate(rolling%correlation(m,m,horizon,0))
   end subroutine allocate_empty_dcc_rolling

end module rmgarch_rolling
