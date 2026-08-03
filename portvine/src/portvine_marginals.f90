! SPDX-License-Identifier: GPL-3.0-only
module portvine_marginals
   use portvine_kinds, only : dp
   use portvine_types, only : marginal_settings_type, asset_marginal_result, &
      marginal_window_result, arma_fit_result, portvine_success, &
      portvine_invalid_input, portvine_fit_failure
   use rugarch, only : garch_spec, garch_fit_result, fit_model, garch_filter, &
      distribution_cdf, absolute_moment, solve_linear_system, &
      model_sgarch, model_gjrgarch, model_igarch, model_egarch, model_aparch
   implicit none
   private
   public :: fit_arma_css, fit_rolling_marginals, get_marginal_point
   public :: marginal_window_number

contains

   subroutine fit_arma_css(y, p, q, fit, max_iterations, tolerance)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: p, q
      type(arma_fit_result), intent(out) :: fit
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: residual(:), old_residual(:), design(:,:), target(:)
      real(dp), allocatable :: xtx(:,:), xty(:), beta(:), old_beta(:)
      real(dp) :: tol, scale
      integer :: n, k, nrow, t, j, iter, maxit, info

      n = size(y)
      k = 1+p+q
      maxit = 60
      if (present(max_iterations)) maxit = max(1,max_iterations)
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = max(tiny(1.0_dp),tolerance)
      allocate(fit%ar(max(0,p)),fit%ma(max(0,q)),fit%residuals(n))
      if (p > 0) fit%ar = 0.0_dp
      if (q > 0) fit%ma = 0.0_dp
      fit%intercept = sum(y)/real(max(1,n),dp)
      fit%residuals = y-fit%intercept
      if (n <= max(p,q)+k+2) then
         fit%status = 2
         return
      end if
      nrow = n-max(p,q)
      allocate(residual(n),old_residual(n),design(nrow,k),target(nrow), &
         xtx(k,k),xty(k),beta(k),old_beta(k))
      residual = fit%residuals
      old_beta = 0.0_dp
      old_beta(1) = fit%intercept
      do iter = 1, maxit
         old_residual = residual
         do t = max(p,q)+1, n
            design(t-max(p,q),1) = 1.0_dp
            do j = 1, p
               design(t-max(p,q),1+j) = y(t-j)
            end do
            do j = 1, q
               design(t-max(p,q),1+p+j) = old_residual(t-j)
            end do
            target(t-max(p,q)) = y(t)
         end do
         xtx = matmul(transpose(design),design)
         scale = max(1.0_dp,maxval(abs([(xtx(j,j),j=1,k)])))
         do j = 1, k
            xtx(j,j) = xtx(j,j)+1.0e-9_dp*scale
         end do
         xty = matmul(transpose(design),target)
         call solve_linear_system(xtx,xty,beta,info)
         if (info /= 0) then
            fit%status = 3
            return
         end if
         beta(2:) = max(-0.995_dp,min(0.995_dp,beta(2:)))
         residual = 0.0_dp
         do t = 1, n
            residual(t) = y(t)-beta(1)
            do j = 1, min(p,t-1)
               residual(t) = residual(t)-beta(1+j)*y(t-j)
            end do
            do j = 1, min(q,t-1)
               residual(t) = residual(t)-beta(1+p+j)*residual(t-j)
            end do
         end do
         if (maxval(abs(beta-old_beta)) <= tol*(1.0_dp+maxval(abs(old_beta)))) exit
         old_beta = beta
      end do
      fit%intercept = beta(1)
      if (p > 0) fit%ar = beta(2:1+p)
      if (q > 0) fit%ma = beta(2+p:1+p+q)
      fit%residuals = residual
      fit%css = sum(residual(max(p,q)+1:)**2)
      fit%iterations = iter
      fit%status = 0
   end subroutine fit_arma_css

   subroutine fit_rolling_marginals(returns, settings, vine_train_size, result, status)
      real(dp), intent(in) :: returns(:,:)
      type(marginal_settings_type), intent(in) :: settings
      integer, intent(in) :: vine_train_size
      type(asset_marginal_result), allocatable, intent(out) :: result(:)
      integer, intent(out), optional :: status
      integer :: n_assets, nobs, nw, a, w, istat
      type(garch_spec) :: base_spec

      n_assets = size(returns,1)
      nobs = size(returns,2)
      if (present(status)) status = portvine_success
      if (n_assets < 1 .or. settings%train_size < 30 .or. &
          settings%train_size >= nobs .or. settings%refit_size < 1 .or. &
          vine_train_size < 3 .or. vine_train_size > settings%train_size .or. &
          .not.allocated(settings%spec) .or. &
          (size(settings%spec) /= 1 .and. size(settings%spec) /= n_assets)) then
         allocate(result(0))
         if (present(status)) status = portvine_invalid_input
         return
      end if
      nw = (nobs-settings%train_size+settings%refit_size-1)/settings%refit_size
      allocate(result(n_assets))
      do a = 1, n_assets
         allocate(result(a)%window(nw))
         if (size(settings%spec) == 1) then
            base_spec = settings%spec(1)
         else
            base_spec = settings%spec(a)
         end if
         result(a)%status = portvine_success
         do w = 1, nw
            call fit_one_window(returns(a,:),base_spec,settings%train_size, &
               settings%refit_size,vine_train_size,w,settings%max_iterations, &
               result(a)%window(w),istat)
            if (istat /= portvine_success) result(a)%status = portvine_fit_failure
         end do
      end do
      if (present(status)) then
         if (any([(result(a)%status /= portvine_success,a=1,n_assets)])) &
            status = portvine_fit_failure
      end if
   end subroutine fit_rolling_marginals

   subroutine fit_one_window(y, base_spec, train_size, refit_size, vine_train_size, &
      window_number, max_iterations, result, status)
      real(dp), intent(in) :: y(:)
      type(garch_spec), intent(in) :: base_spec
      integer, intent(in) :: train_size, refit_size, vine_train_size
      integer, intent(in) :: window_number, max_iterations
      type(marginal_window_result), intent(out) :: result
      integer, intent(out) :: status
      type(arma_fit_result) :: arma
      type(garch_fit_result) :: gfit
      type(garch_spec) :: spec
      real(dp), allocatable :: train(:), full(:), eps(:), sig(:)
      integer :: forecast_start, forecast_end, train_start, train_end, store_start
      integer :: p, q, pa, qm, ntrain, nfull, offset, i
      logical :: valid

      status = portvine_success
      forecast_start = train_size+(window_number-1)*refit_size+1
      forecast_end = min(size(y),train_size+window_number*refit_size)
      train_end = forecast_start-1
      train_start = max(1,train_end-train_size+1)
      store_start = max(train_start,forecast_start-vine_train_size)
      result%start_index = store_start
      result%forecast_start = forecast_start
      result%forecast_end = forecast_end
      if (forecast_start > size(y) .or. train_end-train_start+1 < 30) then
         result%status = portvine_invalid_input
         status = result%status
         return
      end if
      pa = 0
      qm = 0
      p = 0
      q = 0
      if (allocated(base_spec%ar)) pa = size(base_spec%ar)
      if (allocated(base_spec%ma)) qm = size(base_spec%ma)
      if (allocated(base_spec%alpha)) p = size(base_spec%alpha)
      if (allocated(base_spec%beta)) q = size(base_spec%beta)
      allocate(train(train_end-train_start+1))
      train = y(train_start:train_end)
      call fit_arma_css(train,pa,qm,arma,max_iterations=min(80,max_iterations))
      if (arma%status /= 0) then
         result%status = portvine_fit_failure
         status = result%status
         return
      end if
      gfit = fit_model(arma%residuals,base_spec%model,p,q,base_spec%cond_dist, &
         fit_mean=.false.,fit_shape=.true.,fit_skew=.true.,fit_lambda=.true., &
         max_iterations=max_iterations)
      if (gfit%status > 1) then
         result%status = portvine_fit_failure
         status = result%status
         return
      end if
      spec = gfit%spec
      spec%mean = arma%intercept
      if (allocated(spec%ar)) deallocate(spec%ar)
      if (allocated(spec%ma)) deallocate(spec%ma)
      allocate(spec%ar(pa),spec%ma(qm))
      if (pa > 0) spec%ar = arma%ar
      if (qm > 0) spec%ma = arma%ma

      nfull = forecast_end-train_start+1
      ntrain = train_end-train_start+1
      allocate(full(nfull),eps(nfull),sig(nfull))
      full = y(train_start:forecast_end)
      call garch_filter(full(1:ntrain),spec,eps(1:ntrain),sig(1:ntrain),valid)
      if (.not.valid) then
         result%status = portvine_fit_failure
         status = result%status
         return
      end if
      do i = ntrain+1, nfull
         call update_fixed_model(full,spec,eps,sig,i)
      end do
      offset = store_start-train_start
      allocate(result%residual(forecast_end-store_start+1), &
         result%standardized(forecast_end-store_start+1), &
         result%uniform(forecast_end-store_start+1), &
         result%mu(forecast_end-store_start+1), &
         result%sigma(forecast_end-store_start+1))
      result%residual = eps(offset+1:nfull)
      result%sigma = sig(offset+1:nfull)
      result%standardized = result%residual/max(result%sigma,1.0e-12_dp)
      result%mu = full(offset+1:nfull)-result%residual
      do i = 1, size(result%uniform)
         result%uniform(i) = distribution_cdf(result%standardized(i), &
            spec%cond_dist,spec%shape,spec%skew,spec%lambda)
         result%uniform(i) = max(1.0e-10_dp,min(1.0_dp-1.0e-10_dp,result%uniform(i)))
      end do
      result%garch_fit = gfit
      result%combined_spec = spec
      result%status = portvine_success
   end subroutine fit_one_window

   subroutine update_fixed_model(y,spec,residual,sigma,i)
      real(dp), intent(in) :: y(:)
      type(garch_spec), intent(in) :: spec
      real(dp), intent(inout) :: residual(:),sigma(:)
      integer, intent(in) :: i
      integer :: j, pa, qm, p, q
      real(dp) :: h, hp, logh, z, mabs, shock

      pa = 0; qm = 0; p = 0; q = 0
      if (allocated(spec%ar)) pa = size(spec%ar)
      if (allocated(spec%ma)) qm = size(spec%ma)
      if (allocated(spec%alpha)) p = size(spec%alpha)
      if (allocated(spec%beta)) q = size(spec%beta)
      residual(i) = y(i)-spec%mean
      do j = 1, min(pa,i-1)
         residual(i) = residual(i)-spec%ar(j)*(y(i-j)-spec%mean)
      end do
      do j = 1, min(qm,i-1)
         residual(i) = residual(i)-spec%ma(j)*residual(i-j)
      end do
      select case (spec%model)
      case (model_egarch)
         mabs = absolute_moment(1.0_dp,spec%cond_dist,spec%shape,spec%skew,spec%lambda)
         logh = spec%omega
         do j = 1, min(p,i-1)
            z = residual(i-j)/max(sigma(i-j),1.0e-12_dp)
            logh = logh+spec%alpha(j)*(abs(z)-mabs)+spec%gamma(j)*z
         end do
         do j = 1, min(q,i-1)
            logh = logh+spec%beta(j)*log(max(sigma(i-j)**2,1.0e-20_dp))
         end do
         sigma(i) = exp(0.5_dp*max(-100.0_dp,min(100.0_dp,logh)))
      case (model_aparch)
         hp = spec%omega
         do j = 1, min(p,i-1)
            shock = abs(residual(i-j))-spec%gamma(j)*residual(i-j)
            hp = hp+spec%alpha(j)*max(shock,0.0_dp)**spec%delta
         end do
         do j = 1, min(q,i-1)
            hp = hp+spec%beta(j)*sigma(i-j)**spec%delta
         end do
         sigma(i) = max(hp,1.0e-20_dp)**(1.0_dp/spec%delta)
      case default
         h = spec%omega
         do j = 1, min(p,i-1)
            h = h+spec%alpha(j)*residual(i-j)**2
            if (spec%model == model_gjrgarch .and. residual(i-j) < 0.0_dp) &
               h = h+spec%gamma(j)*residual(i-j)**2
         end do
         do j = 1, min(q,i-1)
            h = h+spec%beta(j)*sigma(i-j)**2
         end do
         if (spec%model == model_igarch .and. p+q > 0) h = max(h,1.0e-20_dp)
         sigma(i) = sqrt(max(h,1.0e-20_dp))
      end select
   end subroutine update_fixed_model

   pure integer function marginal_window_number(global_index, train_size, refit_size, nwindow) result(w)
      integer, intent(in) :: global_index, train_size, refit_size, nwindow
      w = (max(global_index,train_size+1)-train_size-1)/refit_size+1
      w = max(1,min(nwindow,w))
   end function marginal_window_number

   subroutine get_marginal_point(asset, window_number, global_index, mu, sigma, uniform, spec, status)
      type(asset_marginal_result), intent(in) :: asset
      integer, intent(in) :: window_number, global_index
      real(dp), intent(out) :: mu, sigma, uniform
      type(garch_spec), intent(out), optional :: spec
      integer, intent(out), optional :: status
      integer :: k
      mu = 0.0_dp
      sigma = 1.0_dp
      uniform = 0.5_dp
      if (present(status)) status = portvine_success
      if (.not.allocated(asset%window) .or. window_number < 1 .or. &
          window_number > size(asset%window)) then
         if (present(status)) status = portvine_invalid_input
         return
      end if
      k = global_index-asset%window(window_number)%start_index+1
      if (k < 1 .or. k > size(asset%window(window_number)%uniform)) then
         if (present(status)) status = portvine_invalid_input
         return
      end if
      mu = asset%window(window_number)%mu(k)
      sigma = asset%window(window_number)%sigma(k)
      uniform = asset%window(window_number)%uniform(k)
      if (present(spec)) spec = asset%window(window_number)%combined_spec
   end subroutine get_marginal_point

end module portvine_marginals
