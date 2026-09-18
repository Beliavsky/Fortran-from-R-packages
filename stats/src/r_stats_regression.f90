! SPDX-License-Identifier: MIT
! SPDX-FileComment: Dense linear-model routines adapted from the R-to-Fortran runtime.
module r_stats_regression
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_distributions, only: r_qt
   use r_kinds, only: dp
   use r_linalg, only: inverse_matrix, least_squares
   use r_stats_design, only: independent_columns
   use r_stats_types, only: lm_fit_t
   implicit none
   private

   public :: lm_aic, lm_coef, lm_confint, lm_cooks_distance, lm_fit, lm_fit_general, &
             lm_predict_general, lm_r_squared_general, predict_lm

   interface lm_fit
      module procedure lm_fit_general
   end interface lm_fit

   interface predict_lm
      module procedure lm_predict_general
   end interface predict_lm

contains

   function lm_fit_general(y, xpred, intercept, weights, rank_tolerance) result(fit)
      !! Fits a weighted dense real linear model with explicit rank-deficiency handling.
      real(dp), intent(in) :: y(:) !! Response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative prior weights with size `n`.
      real(dp), intent(in), optional :: rank_tolerance !! Relative design-rank tolerance.
      type(lm_fit_t) :: fit !! Fitted model and retained fitting data.
      real(dp), allocatable :: design(:, :), design_weighted(:, :), reduced_design(:, :)
      real(dp), allocatable :: reduced_coef(:), reduced_cov(:, :), response_weighted(:), xtx(:, :)
      real(dp) :: rank_tol, sse, sst, weight_sum, ybar
      integer, allocatable :: independent(:)
      integer :: info, j, k, n, n_positive, p

      n = size(y)
      p = size(xpred, 2)
      if (size(xpred, 1) /= n) error stop "lm_fit: response and predictor row counts differ"
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) error stop "lm_fit: invalid weights"
         if (.not. all(ieee_is_finite(weights))) error stop "lm_fit: nonfinite weights"
         fit%weights = weights
         fit%weighted = .true.
      else
         allocate (fit%weights(n), source=1.0_dp)
      end if

      fit%has_intercept = .true.
      if (present(intercept)) fit%has_intercept = intercept
      k = p + merge(1, 0, fit%has_intercept)
      n_positive = count(fit%weights > 0.0_dp)
      if (n_positive == 0) error stop "lm_fit: no positive-weight observations"

      allocate (design(n, k))
      if (fit%has_intercept) then
         design(:, 1) = 1.0_dp
         if (p > 0) design(:, 2:) = xpred
      else if (p > 0) then
         design = xpred
      end if

      design_weighted = design
      response_weighted = y
      do j = 1, n
         design_weighted(j, :) = sqrt(fit%weights(j))*design_weighted(j, :)
         response_weighted(j) = sqrt(fit%weights(j))*response_weighted(j)
      end do
      rank_tol = sqrt(epsilon(1.0_dp))
      if (present(rank_tolerance)) rank_tol = rank_tolerance
      if (rank_tol <= 0.0_dp) error stop "lm_fit: rank tolerance must be positive"
      call independent_columns(design_weighted, rank_tol, independent)
      fit%rank = size(independent)
      if (fit%rank == 0) error stop "lm_fit: zero-rank design matrix"
      reduced_design = design_weighted(:, independent)
      allocate (reduced_coef(fit%rank))
      call least_squares(reduced_design, response_weighted, reduced_coef, info)
      if (info /= 0) error stop "lm_fit: least-squares solution failed"
      allocate (fit%coef(k), source=0.0_dp)
      fit%coef(independent) = reduced_coef
      allocate (fit%aliased(k), source=.true.)
      fit%aliased(independent) = .false.

      xtx = matmul(transpose(reduced_design), reduced_design)
      call inverse_matrix(xtx, reduced_cov, info)
      if (info /= 0) error stop "lm_fit: singular coefficient covariance"
      allocate (fit%cov_unscaled(k, k), source=0.0_dp)
      do j = 1, fit%rank
         fit%cov_unscaled(independent, independent(j)) = reduced_cov(:, j)
      end do

      fit%y = y
      fit%xpred = xpred
      fit%fitted = matmul(design, fit%coef)
      fit%resid = y - fit%fitted
      fit%df = n_positive - fit%rank

      sse = sum(fit%weights*fit%resid**2)
      fit%sigma = sqrt(sse/real(max(1, fit%df), dp))
      if (n > 0) then
         weight_sum = sum(fit%weights)
         if (fit%has_intercept) then
            ybar = sum(fit%weights*y)/weight_sum
            sst = sum(fit%weights*(y - ybar)**2)
         else
            sst = sum(fit%weights*y**2)
         end if
      else
         sst = 0.0_dp
      end if
      if (sst > 0.0_dp) then
         fit%r_squared = 1.0_dp - sse/sst
      else
         fit%r_squared = 0.0_dp
      end if
      if (fit%df > 0 .and. n_positive > 1) then
         fit%adj_r_squared = 1.0_dp - (1.0_dp - fit%r_squared)* &
                             real(n_positive - merge(1, 0, fit%has_intercept), dp)/real(fit%df, dp)
      else
         fit%adj_r_squared = fit%r_squared
      end if
   end function lm_fit_general

   pure function lm_predict_general(fit, xpred) result(yhat)
      !! Predicts responses from a fitted linear model.
      type(lm_fit_t), intent(in) :: fit !! Previously fitted linear model.
      real(dp), intent(in) :: xpred(:, :) !! New predictor rows, excluding an intercept column.
      real(dp), allocatable :: yhat(:) !! Predicted responses with one value per predictor row.
      integer :: p

      p = size(xpred, 2)
      if (fit%has_intercept) then
         if (size(fit%coef) /= p + 1) error stop "predict_lm: predictor count mismatch"
         yhat = fit%coef(1) + matmul(xpred, fit%coef(2:))
      else
         if (size(fit%coef) /= p) error stop "predict_lm: predictor count mismatch"
         yhat = matmul(xpred, fit%coef)
      end if
   end function lm_predict_general

   function lm_coef(y, xpred, intercept) result(coef)
      !! Fits a linear model and returns only its coefficient vector.
      real(dp), intent(in) :: y(:) !! Response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix with shape `(n, p)`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      real(dp), allocatable :: coef(:) !! Estimated coefficient vector.
      type(lm_fit_t) :: fit

      if (present(intercept)) then
         fit = lm_fit_general(y, xpred, intercept)
      else
         fit = lm_fit_general(y, xpred)
      end if
      coef = fit%coef
   end function lm_coef

   function lm_r_squared_general(y, xpred, intercept) result(value)
      !! Fits a linear model and returns its coefficient of determination.
      real(dp), intent(in) :: y(:) !! Response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix with shape `(n, p)`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      real(dp) :: value !! Coefficient of determination.
      type(lm_fit_t) :: fit

      if (present(intercept)) then
         fit = lm_fit_general(y, xpred, intercept)
      else
         fit = lm_fit_general(y, xpred)
      end if
      value = fit%r_squared
   end function lm_r_squared_general

   pure function lm_confint(fit, level) result(intervals)
      !! Computes two-sided coefficient confidence intervals using Student-t quantiles.
      type(lm_fit_t), intent(in) :: fit !! Previously fitted linear model.
      real(dp), intent(in), optional :: level !! Confidence level in `(0, 1)`; defaults to 0.95.
      real(dp), allocatable :: intervals(:, :) !! Lower and upper bounds with shape `(p, 2)`.
      real(dp) :: confidence, critical, standard_error
      integer :: j, p

      confidence = 0.95_dp
      if (present(level)) confidence = level
      if (confidence <= 0.0_dp .or. confidence >= 1.0_dp) then
         error stop "lm_confint: level must lie strictly between zero and one"
      end if
      p = size(fit%coef)
      allocate (intervals(p, 2))
      critical = r_qt(0.5_dp*(1.0_dp + confidence), real(fit%df, dp))
      do j = 1, p
         if (fit%aliased(j)) then
            intervals(j, :) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            standard_error = fit%sigma*sqrt(max(0.0_dp, fit%cov_unscaled(j, j)))
            intervals(j, 1) = fit%coef(j) - critical*standard_error
            intervals(j, 2) = fit%coef(j) + critical*standard_error
         end if
      end do
   end function lm_confint

   pure function lm_cooks_distance(fit) result(distance)
      !! Computes Cook's distance for every observation in a fitted linear model.
      type(lm_fit_t), intent(in) :: fit !! Previously fitted linear model.
      real(dp), allocatable :: distance(:) !! Cook's distances with size `n`.
      real(dp), allocatable :: row(:)
      real(dp) :: denominator, leverage
      integer :: i, k, n, p

      n = size(fit%resid)
      p = fit%rank
      k = size(fit%coef)
      allocate (distance(n), row(k))
      distance = 0.0_dp
      if (p == 0 .or. fit%sigma <= 0.0_dp) return
      do i = 1, n
         if (fit%has_intercept) then
            row(1) = 1.0_dp
            if (p > 1) row(2:) = fit%xpred(i, :)
         else
            row = fit%xpred(i, :)
         end if
         row = sqrt(fit%weights(i))*row
         leverage = dot_product(row, matmul(fit%cov_unscaled, row))
         denominator = real(p, dp)*fit%sigma**2*max(tiny(1.0_dp), (1.0_dp - leverage)**2)
         distance(i) = fit%weights(i)*fit%resid(i)**2*leverage/denominator
      end do
   end function lm_cooks_distance

   pure function lm_aic(fit, penalty) result(value)
      !! Computes the constant-free Gaussian AIC-style score used by stepwise selection.
      type(lm_fit_t), intent(in) :: fit !! Previously fitted linear model.
      real(dp), intent(in), optional :: penalty !! Penalty per coefficient; defaults to two.
      real(dp) :: value !! AIC-style score; lower values indicate preferred models.
      real(dp) :: k, n, rss

      k = 2.0_dp
      if (present(penalty)) k = penalty
      n = real(size(fit%y), dp)
      rss = max(tiny(1.0_dp), sum(fit%resid**2))
      value = n*log(rss/max(1.0_dp, n)) + k*real(fit%rank, dp)
   end function lm_aic

end module r_stats_regression
