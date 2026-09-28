! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
module ftsa_gaev
   use r_kinds, only : dp
   use r_linalg, only : solve_spd, symmetric_eigen
   use ftsa_statistics, only : ftsa_invalid_input, ftsa_ok
   implicit none
   private

   public :: ftsa_gaev_forecast

contains

   subroutine ftsa_gaev_forecast(data, q, d_loc_max, d_logscale_max, density_pred, loc_pred, scale_pred, &
      shape_pred, selected_loc, selected_logscale, status, max_var_lag)
      real(dp), intent(in) :: data(:, :) !! Time-by-grid positive curves used to fit smooth nonstationary GEV models.
      real(dp), intent(in) :: q(:) !! Quantiles/observed values at which forecast GEV densities are evaluated.
      integer, intent(in) :: d_loc_max !! Maximum location spline-basis dimension; values from three through this limit are tested.
      integer, intent(in) :: d_logscale_max !! Maximum log-scale spline dimension; values from three through this limit are tested.
      real(dp), allocatable, intent(out) :: density_pred(:, :) !! `size(q)`-by-grid one-step GEV density forecasts.
      real(dp), allocatable, intent(out) :: loc_pred(:) !! Forecast smooth location parameter on the grid.
      real(dp), allocatable, intent(out) :: scale_pred(:) !! Forecast positive smooth GEV scale on the grid.
      real(dp), intent(out) :: shape_pred !! Forecast constant GEV shape parameter.
      integer, intent(out) :: selected_loc !! Cross-validation-selected location basis dimension.
      integer, intent(out) :: selected_logscale !! Cross-validation-selected log-scale basis dimension.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or failed model fitting.
      integer, intent(in), optional :: max_var_lag !! Maximum lag for coefficient VAR selection; default is three.
      real(dp), allocatable :: basis_loc(:, :), basis_scale(:, :), loc_coef(:, :), logscale_coef(:, :), shape(:)
      real(dp), allocatable :: candidate_loc(:), candidate_scale(:), mode(:)
      real(dp) :: best_loss, loss, candidate_shape
      integer :: k1, k2, local_status, max_lag, n_grid, n_time

      allocate(density_pred(0, 0), loc_pred(0), scale_pred(0))
      selected_loc = 0
      selected_logscale = 0
      shape_pred = 0.0_dp
      n_time = size(data, 1)
      n_grid = size(data, 2)
      max_lag = 3
      if (present(max_var_lag)) max_lag = max_var_lag
      if (n_time < 6 .or. n_grid < 5 .or. size(q) < 1 .or. d_loc_max < 3 .or. d_logscale_max < 3 .or. &
         d_loc_max > n_grid .or. d_logscale_max > n_grid .or. max_lag < 0) then
         status = ftsa_invalid_input
         return
      end if

      best_loss = huge(1.0_dp)
      do k1 = 3, d_loc_max
         do k2 = 3, d_logscale_max
            call gaev_fit_series(data(1:n_time - 1, :), k1, k2, basis_loc, basis_scale, loc_coef, &
               logscale_coef, shape, local_status)
            if (local_status /= ftsa_ok) cycle
            call gaev_forecast_parameters(basis_loc, basis_scale, loc_coef, logscale_coef, shape, max_lag, &
               candidate_loc, candidate_scale, candidate_shape, local_status)
            if (local_status /= ftsa_ok) cycle
            allocate(mode(n_grid))
            call gev_mode_vector(candidate_loc, candidate_scale, candidate_shape, mode)
            loss = symmetric_kld(mode, data(n_time, :))
            deallocate(mode)
            if (loss < best_loss) then
               best_loss = loss
               selected_loc = k1
               selected_logscale = k2
            end if
         end do
      end do
      if (selected_loc < 3 .or. selected_logscale < 3) then
         status = ftsa_invalid_input
         return
      end if

      call gaev_fit_series(data, selected_loc, selected_logscale, basis_loc, basis_scale, loc_coef, &
         logscale_coef, shape, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      call gaev_forecast_parameters(basis_loc, basis_scale, loc_coef, logscale_coef, shape, max_lag, &
         loc_pred, scale_pred, shape_pred, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      deallocate(density_pred)
      allocate(density_pred(size(q), n_grid))
      do k1 = 1, size(q)
         do k2 = 1, n_grid
            density_pred(k1, k2) = gev_pdf(q(k1), loc_pred(k2), scale_pred(k2), shape_pred)
         end do
      end do
      status = ftsa_ok
   end subroutine ftsa_gaev_forecast

   subroutine gaev_fit_series(data, k_loc, k_scale, basis_loc, basis_scale, loc_coef, logscale_coef, shape, status)
      real(dp), intent(in) :: data(:, :) !! Time-by-grid observations to fit independently with smooth GEV parameter curves.
      integer, intent(in) :: k_loc !! Number of B-spline basis columns for the location function.
      integer, intent(in) :: k_scale !! Number of B-spline basis columns for the log-scale function.
      real(dp), allocatable, intent(out) :: basis_loc(:, :) !! Grid-by-location-basis design matrix.
      real(dp), allocatable, intent(out) :: basis_scale(:, :) !! Grid-by-log-scale-basis design matrix.
      real(dp), allocatable, intent(out) :: loc_coef(:, :) !! Time-by-`k_loc` fitted location coefficients.
      real(dp), allocatable, intent(out) :: logscale_coef(:, :) !! Time-by-`k_scale` fitted log-scale coefficients.
      real(dp), allocatable, intent(out) :: shape(:) !! Time-varying fitted constant GEV shape values.
      integer, intent(out) :: status !! Zero when every time slice is fitted successfully.
      real(dp), allocatable :: beta_loc(:), beta_scale(:)
      real(dp) :: shape_value
      integer :: local_status, t

      allocate(basis_loc(0, 0), basis_scale(0, 0), loc_coef(0, 0), logscale_coef(0, 0), shape(0))
      call bspline_basis(size(data, 2), k_loc, basis_loc, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      call bspline_basis(size(data, 2), k_scale, basis_scale, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      deallocate(loc_coef, logscale_coef, shape)
      allocate(loc_coef(size(data, 1), k_loc), logscale_coef(size(data, 1), k_scale), shape(size(data, 1)))
      do t = 1, size(data, 1)
         call fit_smooth_gev(data(t, :), basis_loc, basis_scale, beta_loc, beta_scale, shape_value, local_status)
         if (local_status /= ftsa_ok) then
            status = local_status
            return
         end if
         loc_coef(t, :) = beta_loc
         logscale_coef(t, :) = beta_scale
         shape(t) = shape_value
      end do
      status = ftsa_ok
   end subroutine gaev_fit_series

   subroutine fit_smooth_gev(y, basis_loc, basis_scale, beta_loc, beta_scale, shape_value, status)
      real(dp), intent(in) :: y(:) !! One functional observation across the grid for smooth GEV fitting.
      real(dp), intent(in) :: basis_loc(:, :) !! Grid-by-location spline basis.
      real(dp), intent(in) :: basis_scale(:, :) !! Grid-by-log-scale spline basis.
      real(dp), allocatable, intent(out) :: beta_loc(:) !! Fitted location basis coefficients.
      real(dp), allocatable, intent(out) :: beta_scale(:) !! Fitted log-scale basis coefficients.
      real(dp), intent(out) :: shape_value !! Fitted bounded GEV shape parameter.
      integer, intent(out) :: status !! Zero if the penalized likelihood optimization returns finite parameters.
      real(dp), allocatable :: beta_init(:), gradient(:), hessian_inv(:, :), old_gradient(:), p(:), theta(:), trial(:)
      real(dp), allocatable :: s(:), ydiff(:)
      real(dp) :: alpha, denom, f, ftrial, range_y, scale0, slope
      integer :: d, i, info, iter, k_loc, k_scale

      allocate(beta_loc(0), beta_scale(0))
      k_loc = size(basis_loc, 2)
      k_scale = size(basis_scale, 2)
      if (size(y) /= size(basis_loc, 1) .or. size(y) /= size(basis_scale, 1) .or. k_loc < 1 .or. k_scale < 1) then
         status = ftsa_invalid_input
         return
      end if
      call ridge_least_squares_vector(basis_loc, y, beta_init, info)
      if (info /= 0) then
         status = ftsa_invalid_input
         return
      end if
      range_y = maxval(y) - minval(y)
      scale0 = max(sqrt(sum((y - sum(y) / real(size(y), dp))**2) / real(max(1, size(y) - 1), dp)), &
         max(1.0e-3_dp, 0.05_dp * max(1.0_dp, range_y)))
      d = k_loc + k_scale + 1
      allocate(theta(d), gradient(d), old_gradient(d), hessian_inv(d, d), p(d), trial(d), s(d), ydiff(d))
      theta(1:k_loc) = beta_init
      theta(k_loc + 1:k_loc + k_scale) = log(scale0)
      theta(d) = 0.0_dp
      hessian_inv = 0.0_dp
      do i = 1, d
         hessian_inv(i, i) = 1.0_dp
      end do
      f = gev_objective(y, basis_loc, basis_scale, theta)
      call gev_gradient(y, basis_loc, basis_scale, theta, gradient)
      do iter = 1, 80
         if (maxval(abs(gradient)) < 2.0e-5_dp) exit
         p = -matmul(hessian_inv, gradient)
         slope = dot_product(gradient, p)
         if (slope >= -tiny(1.0_dp)) then
            p = -gradient
            slope = -dot_product(gradient, gradient)
            hessian_inv = 0.0_dp
            do i = 1, d
               hessian_inv(i, i) = 1.0_dp
            end do
         end if
         alpha = 1.0_dp
         do i = 1, 24
            trial = theta + alpha * p
            ftrial = gev_objective(y, basis_loc, basis_scale, trial)
            if (ftrial < huge(1.0_dp) / 100.0_dp .and. ftrial <= f + 1.0e-4_dp * alpha * slope) exit
            alpha = 0.5_dp * alpha
         end do
         if (alpha < 1.0e-7_dp) exit
         s = trial - theta
         old_gradient = gradient
         theta = trial
         f = ftrial
         call gev_gradient(y, basis_loc, basis_scale, theta, gradient)
         ydiff = gradient - old_gradient
         denom = dot_product(ydiff, s)
         if (denom > 1.0e-10_dp * sqrt(max(tiny(1.0_dp), dot_product(s, s) * dot_product(ydiff, ydiff)))) then
            call bfgs_inverse_update(hessian_inv, s, ydiff, denom)
         else
            hessian_inv = 0.0_dp
            do i = 1, d
               hessian_inv(i, i) = 1.0_dp
            end do
         end if
      end do
      if (f >= huge(1.0_dp) / 100.0_dp) then
         status = ftsa_invalid_input
         return
      end if
      deallocate(beta_loc, beta_scale)
      allocate(beta_loc(k_loc), beta_scale(k_scale))
      beta_loc = theta(1:k_loc)
      beta_scale = theta(k_loc + 1:k_loc + k_scale)
      shape_value = 0.45_dp * tanh(theta(d))
      status = ftsa_ok
   end subroutine fit_smooth_gev

   real(dp) function gev_objective(y, basis_loc, basis_scale, theta) result(value)
      real(dp), intent(in) :: y(:) !! Grid observations entering the smooth GEV negative log likelihood.
      real(dp), intent(in) :: basis_loc(:, :) !! Location spline basis.
      real(dp), intent(in) :: basis_scale(:, :) !! Log-scale spline basis.
      real(dp), intent(in) :: theta(:) !! Location coefficients, log-scale coefficients, and unconstrained shape coordinate.
      real(dp), allocatable :: eta_scale(:), location(:)
      real(dp) :: logpdf, penalty, shape, scale, t, z
      integer :: i, k_loc, k_scale

      k_loc = size(basis_loc, 2)
      k_scale = size(basis_scale, 2)
      if (size(theta) /= k_loc + k_scale + 1) then
         value = huge(1.0_dp)
         return
      end if
      if (any(abs(theta(k_loc + 1:k_loc + k_scale)) > 12.0_dp) .or. any(abs(theta(1:k_loc)) > 1.0e6_dp)) then
         value = huge(1.0_dp)
         return
      end if
      location = matmul(basis_loc, theta(1:k_loc))
      eta_scale = matmul(basis_scale, theta(k_loc + 1:k_loc + k_scale))
      shape = 0.45_dp * tanh(theta(size(theta)))
      value = 0.0_dp
      do i = 1, size(y)
         scale = exp(eta_scale(i))
         z = (y(i) - location(i)) / scale
         if (abs(shape) < 1.0e-6_dp) then
            logpdf = -eta_scale(i) - z - safe_exp(-z)
         else
            t = 1.0_dp + shape * z
            if (t <= 1.0e-10_dp) then
               value = huge(1.0_dp)
               return
            end if
            logpdf = -eta_scale(i) - (1.0_dp + 1.0_dp / shape) * log(t) - &
               safe_exp((-1.0_dp / shape) * log(t))
         end if
         value = value - logpdf
      end do
      penalty = 0.0_dp
      if (k_loc >= 3) penalty = penalty + sum((theta(3:k_loc) - 2.0_dp * theta(2:k_loc - 1) + theta(1:k_loc - 2))**2)
      if (k_scale >= 3) penalty = penalty + &
         sum((theta(k_loc + 3:k_loc + k_scale) - 2.0_dp * theta(k_loc + 2:k_loc + k_scale - 1) + &
         theta(k_loc + 1:k_loc + k_scale - 2))**2)
      value = value + 1.0e-3_dp * penalty
   end function gev_objective

   subroutine gev_gradient(y, basis_loc, basis_scale, theta, gradient)
      real(dp), intent(in) :: y(:) !! Grid observations used by the finite-difference likelihood gradient.
      real(dp), intent(in) :: basis_loc(:, :) !! Location spline basis used in the objective.
      real(dp), intent(in) :: basis_scale(:, :) !! Log-scale spline basis used in the objective.
      real(dp), intent(in) :: theta(:) !! Current unconstrained smooth-GEV parameter vector.
      real(dp), intent(out) :: gradient(size(theta)) !! Central finite-difference gradient of the penalized negative log likelihood.
      real(dp) :: h, left(size(theta)), right(size(theta))
      integer :: i

      do i = 1, size(theta)
         h = 1.0e-5_dp * max(1.0_dp, abs(theta(i)))
         left = theta
         right = theta
         left(i) = left(i) - h
         right(i) = right(i) + h
         gradient(i) = (gev_objective(y, basis_loc, basis_scale, right) - &
            gev_objective(y, basis_loc, basis_scale, left)) / (2.0_dp * h)
      end do
   end subroutine gev_gradient

   subroutine bfgs_inverse_update(hessian_inv, s, ydiff, denom)
      real(dp), intent(inout) :: hessian_inv(:, :) !! Dense inverse-Hessian approximation updated in place.
      real(dp), intent(in) :: s(:) !! Accepted parameter step.
      real(dp), intent(in) :: ydiff(:) !! Change in gradient across the accepted step.
      real(dp), intent(in) :: denom !! Positive curvature product `dot_product(ydiff,s)`.
      real(dp) :: hy(size(s)), yhy, factor
      integer :: i, j

      hy = matmul(hessian_inv, ydiff)
      yhy = dot_product(ydiff, hy)
      factor = (1.0_dp + yhy / denom) / denom
      do j = 1, size(s)
         do i = 1, size(s)
            hessian_inv(i, j) = hessian_inv(i, j) + factor * s(i) * s(j) - &
               (s(i) * hy(j) + hy(i) * s(j)) / denom
         end do
      end do
   end subroutine bfgs_inverse_update

   subroutine gaev_forecast_parameters(basis_loc, basis_scale, loc_coef, logscale_coef, shape, max_lag, &
      loc_pred, scale_pred, shape_pred, status)
      real(dp), intent(in) :: basis_loc(:, :) !! Grid-by-location spline basis used for reconstruction.
      real(dp), intent(in) :: basis_scale(:, :) !! Grid-by-log-scale spline basis used for reconstruction.
      real(dp), intent(in) :: loc_coef(:, :) !! Time-by-location coefficient history.
      real(dp), intent(in) :: logscale_coef(:, :) !! Time-by-log-scale coefficient history.
      real(dp), intent(in) :: shape(:) !! Historical GEV shape estimates.
      integer, intent(in) :: max_lag !! Maximum VAR lag considered by AIC; zero permits trend-only forecasting.
      real(dp), allocatable, intent(out) :: loc_pred(:) !! One-step forecast location curve.
      real(dp), allocatable, intent(out) :: scale_pred(:) !! One-step forecast positive scale curve.
      real(dp), intent(out) :: shape_pred !! One-step forecast GEV shape, bounded to a stable interval.
      integer, intent(out) :: status !! Zero on success; nonzero when no coefficient VAR/trend fit is feasible.
      real(dp), allocatable :: allpara(:, :), coefficients(:, :), residuals(:, :), one(:)
      integer :: k_loc, k_scale, lag, local_status, n

      allocate(loc_pred(0), scale_pred(0))
      k_loc = size(loc_coef, 2)
      k_scale = size(logscale_coef, 2)
      n = size(loc_coef, 1)
      if (size(logscale_coef, 1) /= n .or. size(shape) /= n .or. n < 3) then
         status = ftsa_invalid_input
         return
      end if
      allocate(allpara(n, k_loc + k_scale + 1))
      allpara(:, 1:k_loc) = loc_coef
      allpara(:, k_loc + 1:k_loc + k_scale) = logscale_coef
      allpara(:, size(allpara, 2)) = shape
      call select_var_aic_trend(allpara, max_lag, lag, coefficients, residuals, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      call forecast_var_trend(allpara, coefficients, lag, one, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      deallocate(loc_pred, scale_pred)
      allocate(loc_pred(size(basis_loc, 1)), scale_pred(size(basis_scale, 1)))
      loc_pred = matmul(basis_loc, one(1:k_loc))
      scale_pred = exp(matmul(basis_scale, one(k_loc + 1:k_loc + k_scale)))
      shape_pred = min(0.45_dp, max(-0.45_dp, one(size(one))))
      status = ftsa_ok
   end subroutine gaev_forecast_parameters

   subroutine select_var_aic_trend(series, max_lag, selected_lag, coefficients, residuals, status)
      real(dp), intent(in) :: series(:, :) !! Time-by-variable parameter history used in a VAR with intercept and trend.
      integer, intent(in) :: max_lag !! Largest autoregressive lag considered; zero allows a trend-only model.
      integer, intent(out) :: selected_lag !! Lag minimizing the regularized Gaussian AIC.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Selected regression coefficient matrix.
      real(dp), allocatable, intent(out) :: residuals(:, :) !! Residual matrix from the selected parameter model.
      integer, intent(out) :: status !! Zero when at least one trend/VAR model is estimable.
      real(dp), allocatable :: candidate_coef(:, :), candidate_resid(:, :), covariance(:, :), eval(:), evec(:, :)
      real(dp) :: aic, best, floor_value, logdet
      integer :: d, info, local_status, nobs, npar, p

      allocate(coefficients(0, 0), residuals(0, 0))
      selected_lag = -1
      d = size(series, 2)
      best = huge(1.0_dp)
      do p = 0, max_lag
         call fit_var_trend(series, p, candidate_coef, candidate_resid, local_status)
         if (local_status /= ftsa_ok .or. size(candidate_resid, 1) < 2) cycle
         allocate(covariance(d, d))
         covariance = matmul(transpose(candidate_resid), candidate_resid) / real(size(candidate_resid, 1), dp)
         call symmetric_eigen(0.5_dp * (covariance + transpose(covariance)), eval, evec, info)
         deallocate(covariance)
         if (info /= 0) cycle
         floor_value = 1.0e-10_dp * max(1.0_dp, maxval(abs(eval)))
         logdet = sum(log(max(eval, floor_value)))
         nobs = size(candidate_resid, 1)
         npar = d * (2 + d * p)
         aic = logdet + 2.0_dp * real(npar, dp) / real(nobs, dp)
         if (aic < best) then
            best = aic
            selected_lag = p
            if (allocated(coefficients)) deallocate(coefficients)
            if (allocated(residuals)) deallocate(residuals)
            allocate(coefficients(size(candidate_coef, 1), size(candidate_coef, 2)))
            allocate(residuals(size(candidate_resid, 1), size(candidate_resid, 2)))
            coefficients = candidate_coef
            residuals = candidate_resid
         end if
      end do
      if (selected_lag < 0) then
         status = ftsa_invalid_input
      else
         status = ftsa_ok
      end if
   end subroutine select_var_aic_trend

   subroutine fit_var_trend(series, lag_order, coefficients, residuals, status)
      real(dp), intent(in) :: series(:, :) !! Time-by-variable series fit by multivariate regression with trend and lags.
      integer, intent(in) :: lag_order !! Nonnegative autoregressive lag order.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! `(2+d*lag_order)`-by-`d` fitted coefficient matrix.
      real(dp), allocatable, intent(out) :: residuals(:, :) !! Regression residual matrix after dropping initial lag observations.
      integer, intent(out) :: status !! Zero when the least-squares system is estimable.
      real(dp), allocatable :: design(:, :), response(:, :)
      integer :: d, i, info, lag, n, npar, row, col

      allocate(coefficients(0, 0), residuals(0, 0))
      n = size(series, 1)
      d = size(series, 2)
      npar = 2 + d * lag_order
      if (n - lag_order < max(3, npar) .or. d < 1) then
         status = ftsa_invalid_input
         return
      end if
      allocate(design(n - lag_order, npar), response(n - lag_order, d))
      do row = 1, n - lag_order
         i = lag_order + row
         design(row, 1) = 1.0_dp
         design(row, 2) = real(i, dp)
         col = 2
         do lag = 1, lag_order
            design(row, col + 1:col + d) = series(i - lag, :)
            col = col + d
         end do
         response(row, :) = series(i, :)
      end do
      call ridge_least_squares_matrix(design, response, coefficients, info)
      if (info /= 0) then
         status = ftsa_invalid_input
      else
         deallocate(residuals)
         allocate(residuals(size(response, 1), d))
         residuals = response - matmul(design, coefficients)
         status = ftsa_ok
      end if
   end subroutine fit_var_trend

   subroutine forecast_var_trend(series, coefficients, lag_order, forecast, status)
      real(dp), intent(in) :: series(:, :) !! Historical time-by-variable parameter series.
      real(dp), intent(in) :: coefficients(:, :) !! Coefficients from `fit_var_trend` at the requested lag.
      integer, intent(in) :: lag_order !! Autoregressive lag order represented in `coefficients`.
      real(dp), allocatable, intent(out) :: forecast(:) !! One-step multivariate forecast.
      integer, intent(out) :: status !! Zero when dimensions are consistent.
      real(dp), allocatable :: row(:)
      integer :: col, d, lag, n

      allocate(forecast(0))
      n = size(series, 1)
      d = size(series, 2)
      if (size(coefficients, 1) /= 2 + d * lag_order .or. size(coefficients, 2) /= d) then
         status = ftsa_invalid_input
         return
      end if
      allocate(row(size(coefficients, 1)))
      row(1) = 1.0_dp
      row(2) = real(n + 1, dp)
      col = 2
      do lag = 1, lag_order
         row(col + 1:col + d) = series(n + 1 - lag, :)
         col = col + d
      end do
      deallocate(forecast)
      allocate(forecast(d))
      forecast = matmul(row, coefficients)
      status = ftsa_ok
   end subroutine forecast_var_trend

   subroutine ridge_least_squares_vector(design, response, coefficients, info)
      real(dp), intent(in) :: design(:, :) !! Regression design matrix with observations in rows.
      real(dp), intent(in) :: response(:) !! Response vector with one value per design row.
      real(dp), allocatable, intent(out) :: coefficients(:) !! Ridge-stabilized least-squares coefficient vector.
      integer, intent(out) :: info !! Zero on success; nonzero when dimensions are invalid or the SPD solve fails.
      real(dp), allocatable :: gram(:, :), rhs(:, :), solution(:, :)
      real(dp) :: diagonal_scale, ridge
      integer :: j, npar

      allocate(coefficients(0))
      npar = size(design, 2)
      if (size(design, 1) /= size(response) .or. npar < 1 .or. size(design, 1) < 1) then
         info = 1
         return
      end if
      allocate(gram(npar, npar), rhs(npar, 1))
      gram = matmul(transpose(design), design)
      rhs(:, 1) = matmul(transpose(design), response)
      diagonal_scale = sum(abs([(gram(j, j), j=1, npar)])) / real(npar, dp)
      ridge = sqrt(epsilon(1.0_dp)) * max(1.0_dp, diagonal_scale)
      do j = 1, npar
         gram(j, j) = gram(j, j) + ridge
      end do
      allocate(solution(npar, 1))
      call solve_spd(gram, rhs, solution, info)
      if (info == 0) then
         deallocate(coefficients)
         allocate(coefficients(npar))
         coefficients = solution(:, 1)
      end if
   end subroutine ridge_least_squares_vector

   subroutine ridge_least_squares_matrix(design, response, coefficients, info)
      real(dp), intent(in) :: design(:, :) !! Regression design matrix with observations in rows.
      real(dp), intent(in) :: response(:, :) !! Matrix of simultaneous response series.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Ridge-stabilized least-squares coefficient matrix.
      integer, intent(out) :: info !! Zero on success; nonzero when dimensions are invalid or the SPD solve fails.
      real(dp), allocatable :: gram(:, :), rhs(:, :)
      real(dp) :: diagonal_scale, ridge
      integer :: j, npar

      allocate(coefficients(0, 0))
      npar = size(design, 2)
      if (size(design, 1) /= size(response, 1) .or. npar < 1 .or. size(response, 2) < 1) then
         info = 1
         return
      end if
      allocate(gram(npar, npar), rhs(npar, size(response, 2)))
      gram = matmul(transpose(design), design)
      rhs = matmul(transpose(design), response)
      diagonal_scale = sum(abs([(gram(j, j), j=1, npar)])) / real(npar, dp)
      ridge = sqrt(epsilon(1.0_dp)) * max(1.0_dp, diagonal_scale)
      do j = 1, npar
         gram(j, j) = gram(j, j) + ridge
      end do
      deallocate(coefficients)
      allocate(coefficients(npar, size(response, 2)))
      call solve_spd(gram, rhs, coefficients, info)
   end subroutine ridge_least_squares_matrix

   pure real(dp) function safe_exp(x) result(value)
      real(dp), intent(in) :: x !! Exponent to evaluate while avoiding overflow and subnormal underflow.

      if (x <= -700.0_dp) then
         value = 0.0_dp
      else if (x >= 700.0_dp) then
         value = exp(700.0_dp)
      else
         value = exp(x)
      end if
   end function safe_exp

   subroutine bspline_basis(n_grid, n_basis, basis, status)
      integer, intent(in) :: n_grid !! Number of equally spaced grid points on `[0,1]`.
      integer, intent(in) :: n_basis !! Number of open-uniform B-spline basis functions, at least three.
      real(dp), allocatable, intent(out) :: basis(:, :) !! Grid-by-basis open-uniform spline design matrix.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid basis/grid dimensions.
      real(dp), allocatable :: current(:), next(:), knots(:)
      real(dp) :: denom1, denom2, x
      integer :: degree, i, j, level, n_knots

      allocate(basis(0, 0))
      if (n_grid < 2 .or. n_basis < 3) then
         status = ftsa_invalid_input
         return
      end if
      degree = min(3, n_basis - 1)
      n_knots = n_basis + degree + 1
      allocate(knots(n_knots))
      knots(1:degree + 1) = 0.0_dp
      knots(n_knots - degree:n_knots) = 1.0_dp
      if (n_basis - degree - 1 > 0) then
         do j = 1, n_basis - degree - 1
            knots(degree + 1 + j) = real(j, dp) / real(n_basis - degree, dp)
         end do
      end if
      deallocate(basis)
      allocate(basis(n_grid, n_basis))
      allocate(current(n_basis + degree), next(n_basis + degree))
      do i = 1, n_grid
         x = real(i - 1, dp) / real(n_grid - 1, dp)
         current = 0.0_dp
         if (i == n_grid) then
            current(n_basis) = 1.0_dp
         else
            do j = 1, n_basis + degree
               if (knots(j) <= x .and. x < knots(j + 1)) current(j) = 1.0_dp
            end do
         end if
         do level = 1, degree
            next = 0.0_dp
            do j = 1, n_basis + degree - level
               denom1 = knots(j + level) - knots(j)
               denom2 = knots(j + level + 1) - knots(j + 1)
               if (denom1 > 0.0_dp) next(j) = next(j) + (x - knots(j)) / denom1 * current(j)
               if (denom2 > 0.0_dp) next(j) = next(j) + (knots(j + level + 1) - x) / denom2 * current(j + 1)
            end do
            current = next
         end do
         basis(i, :) = current(1:n_basis)
         if (i == n_grid) then
            basis(i, :) = 0.0_dp
            basis(i, n_basis) = 1.0_dp
         end if
      end do
      status = ftsa_ok
   end subroutine bspline_basis

   pure real(dp) function gev_pdf(x, location, scale, shape) result(pdf)
      real(dp), intent(in) :: x !! Evaluation point for the generalized extreme-value density.
      real(dp), intent(in) :: location !! GEV location parameter.
      real(dp), intent(in) :: scale !! Positive GEV scale parameter.
      real(dp), intent(in) :: shape !! GEV shape parameter; zero denotes the Gumbel limit.
      real(dp) :: t, z

      pdf = 0.0_dp
      if (scale <= 0.0_dp) return
      z = (x - location) / scale
      if (abs(shape) < 1.0e-7_dp) then
         pdf = safe_exp(-z - safe_exp(-z)) / scale
      else
         t = 1.0_dp + shape * z
         if (t <= 0.0_dp) return
         pdf = safe_exp(-(1.0_dp + 1.0_dp / shape) * log(t) - &
            safe_exp((-1.0_dp / shape) * log(t))) / scale
      end if
   end function gev_pdf

   pure subroutine gev_mode_vector(location, scale, shape, mode)
      real(dp), intent(in) :: location(:) !! Gridwise GEV location parameters.
      real(dp), intent(in) :: scale(:) !! Gridwise positive GEV scale parameters.
      real(dp), intent(in) :: shape !! Constant GEV shape parameter.
      real(dp), intent(out) :: mode(size(location)) !! Gridwise modal values used by upstream cross-validation.
      real(dp) :: factor

      if (abs(shape) < 1.0e-7_dp) then
         mode = location
      else
         factor = (sign(1.0_dp, 1.0_dp + shape) * abs(1.0_dp + shape)**(-shape) - 1.0_dp) / shape
         mode = location + scale * factor
      end if
   end subroutine gev_mode_vector

   pure real(dp) function symmetric_kld(px_input, py_input) result(value)
      real(dp), intent(in) :: px_input(:) !! First vector treated as probabilities or log-probabilities as in LaplacesDemon `KLD`.
      real(dp), intent(in) :: py_input(:) !! Second vector treated as probabilities or log-probabilities as in LaplacesDemon `KLD`.
      real(dp) :: px(size(px_input)), py(size(py_input))

      if (size(px_input) /= size(py_input) .or. size(px_input) < 1) then
         value = huge(1.0_dp)
         return
      end if
      px = px_input
      py = py_input
      if (any(px <= 0.0_dp)) px = exp(max(-700.0_dp, min(700.0_dp, px)))
      if (any(py <= 0.0_dp)) py = exp(max(-700.0_dp, min(700.0_dp, py)))
      px = max(px, tiny(1.0_dp))
      py = max(py, tiny(1.0_dp))
      px = px / sum(px)
      py = py / sum(py)
      value = 0.5_dp * (sum(px * log(px / py)) + sum(py * log(py / px)))
   end function symmetric_kld

end module ftsa_gaev
