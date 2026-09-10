! SPDX-License-Identifier: GPL-2.0-only
module marss_residuals_full
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_kf_result, marss_hatyt_result, marss_residual_result
   use marss_analysis, only : marss_hatyt
   use marss_parameters, only : marss_b_at, marss_u_at, marss_q_at, marss_z_at, marss_a_at, marss_r_at
   use marss_covariance, only : psd_inverse, psd_inverse_logdet, psd_cholesky_inverse
   use marss_covariance, only : psd_cholesky_standardize
   implicit none
   private
   public :: marss_residuals_smoothed
   public :: marss_residuals_filtered
   public :: marss_residuals_one_step
   public :: marss_residuals_harvey

contains

   pure subroutine marss_residuals_smoothed(model, kf, result, normalize)
      type(marss_model), intent(in) :: model !! Fitted MARSS model whose residuals are conditioned on the complete data record.
      type(marss_kf_result), intent(in) :: kf !! Filter/smoother output including smoothed lag-one covariances.
      type(marss_residual_result), intent(out) :: result !! Full model/state residual hierarchy corresponding to MARSS tT residuals.
      logical, intent(in), optional :: normalize !! If true, normalize residuals using effective R and Q covariances.
      type(marss_hatyt_result) :: hat
      real(dp) :: at(size(model%a))
      real(dp) :: bnext(size(model%b, 1), size(model%b, 2))
      real(dp) :: cov_et(size(model%z, 1), size(model%b, 1))
      real(dp) :: cov_yx(size(model%z, 1), size(model%b, 1))
      real(dp) :: cov_yxnext(size(model%z, 1), size(model%b, 1))
      real(dp) :: model_var(size(model%z, 1), size(model%z, 1))
      real(dp) :: qnext(size(model%q, 1), size(model%q, 2))
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      real(dp) :: state_var(size(model%b, 1), size(model%b, 1))
      real(dp) :: unext(size(model%u))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      real(dp) :: nan_value
      logical :: do_normalize
      logical :: missing(size(model%z, 1))
      integer :: m
      integer :: n
      integer :: t
      integer :: tt

      call initialize_full_result(model, result)
      call marss_hatyt(model, kf, hat)
      n = size(model%y, 1)
      m = size(model%b, 1)
      tt = size(model%y, 2)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      do_normalize = .false.
      if (present(normalize)) do_normalize = normalize

      do t = 1, tt
         zt = marss_z_at(model, t)
         at = marss_a_at(model, t)
         rt = marss_r_at(model, t)
         missing = ieee_is_nan(model%y(:, t))
         result%residuals(1:n, t) = hat%yt(:, t) - matmul(zt, kf%x_smooth(:, t)) - at
         result%expected_observed_residuals(:, t) = result%residuals(1:n, t)
         result%var_observed_residuals(:, :, t) = hat%ot(:, :, t) - outer_product(hat%yt(:, t), hat%yt(:, t))
         cov_yx = hat%yxt(:, :, t) - outer_product(hat%yt(:, t), kf%x_smooth(:, t))
         model_var = rt - matmul(matmul(zt, kf%p_smooth(:, :, t)), transpose(zt)) + &
            matmul(cov_yx, transpose(zt)) + matmul(zt, transpose(cov_yx))
         model_var = 0.5_dp * (model_var + transpose(model_var))
         result%var_residuals(1:n, 1:n, t) = model_var

         if (t < tt) then
            bnext = marss_b_at(model, t + 1)
            unext = marss_u_at(model, t + 1)
            qnext = marss_q_at(model, t + 1)
            result%residuals(n + 1:n + m, t) = kf%x_smooth(:, t + 1) - &
               matmul(bnext, kf%x_smooth(:, t)) - unext
            cov_yxnext = hat%yxtp(:, :, t) - outer_product(hat%yt(:, t), kf%x_smooth(:, t + 1))
            cov_et = matmul(zt, transpose(kf%p_lag(:, :, t + 1))) - &
               matmul(matmul(zt, kf%p_smooth(:, :, t)), transpose(bnext)) - cov_yxnext + &
               matmul(cov_yx, transpose(bnext))
            state_var = qnext - kf%p_smooth(:, :, t + 1) - &
               matmul(matmul(bnext, kf%p_smooth(:, :, t)), transpose(bnext)) + &
               matmul(kf%p_lag(:, :, t + 1), transpose(bnext)) + &
               matmul(bnext, transpose(kf%p_lag(:, :, t + 1)))
            state_var = 0.5_dp * (state_var + transpose(state_var))
            result%var_residuals(1:n, n + 1:n + m, t) = cov_et
            result%var_residuals(n + 1:n + m, 1:n, t) = transpose(cov_et)
            result%var_residuals(n + 1:n + m, n + 1:n + m, t) = state_var
            if (do_normalize) call normalize_residual(model, t, result%residuals(:, t), result%var_residuals(:, :, t), .true.)
            call standardize_residual(result%residuals(:, t), result%var_residuals(:, :, t), n, missing, .true., &
               result%std_residuals(:, t), result%marginal_residuals(:, t), &
               result%block_cholesky_residuals(:, t), result%info)
            if (result%info /= 0) return
         else
            result%residuals(n + 1:n + m, t) = nan_value
            result%var_residuals(1:n, n + 1:n + m, t) = nan_value
            result%var_residuals(n + 1:n + m, 1:n, t) = nan_value
            result%var_residuals(n + 1:n + m, n + 1:n + m, t) = nan_value
            if (do_normalize) call normalize_observation(model, t, result%residuals(1:n, t), &
               result%var_residuals(1:n, 1:n, t))
            call standardize_residual(result%residuals(1:n, t), result%var_residuals(1:n, 1:n, t), n, missing, .false., &
               result%std_residuals(1:n, t), result%marginal_residuals(1:n, t), &
               result%block_cholesky_residuals(1:n, t), result%info)
            if (result%info /= 0) return
            result%std_residuals(:, t) = nan_value
            result%marginal_residuals(n + 1:n + m, t) = nan_value
            result%block_cholesky_residuals(n + 1:n + m, t) = nan_value
         end if
         where (missing)
            result%residuals(1:n, t) = nan_value
         end where
      end do
      call finalize_full_result(model, kf, result)
   end subroutine marss_residuals_smoothed

   pure subroutine marss_residuals_filtered(model, kf, result, normalize)
      type(marss_model), intent(in) :: model !! Fitted MARSS model whose residuals are conditioned on observations through time t.
      type(marss_kf_result), intent(in) :: kf !! Filter/smoother output; filtered moments are used by this residual variant.
      type(marss_residual_result), intent(out) :: result !! Full output with observation residuals and undefined state residuals.
      logical, intent(in), optional :: normalize !! If true, pre-normalize observation residuals by the effective R square root.
      real(dp) :: at(size(model%a))
      real(dp) :: cross_yx(size(model%z, 1), size(model%b, 1))
      real(dp) :: mean_y(size(model%z, 1))
      real(dp) :: model_var(size(model%z, 1), size(model%z, 1))
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      real(dp) :: second_y(size(model%z, 1), size(model%z, 1))
      real(dp) :: st(size(model%z, 1), size(model%b, 1))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      real(dp) :: nan_value
      logical :: do_normalize
      logical :: missing(size(model%z, 1))
      integer :: m
      integer :: n
      integer :: t
      integer :: tt

      call initialize_full_result(model, result)
      n = size(model%y, 1)
      m = size(model%b, 1)
      tt = size(model%y, 2)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      do_normalize = .false.
      if (present(normalize)) do_normalize = normalize
      do t = 1, tt
         zt = marss_z_at(model, t)
         at = marss_a_at(model, t)
         rt = marss_r_at(model, t)
         missing = ieee_is_nan(model%y(:, t))
         call current_observation_moments(model, t, kf%x_filt(:, t), kf%p_filt(:, :, t), mean_y, second_y, cross_yx)
         result%residuals(1:n, t) = mean_y - matmul(zt, kf%x_filt(:, t)) - at
         result%expected_observed_residuals(:, t) = result%residuals(1:n, t)
         result%var_observed_residuals(:, :, t) = second_y - outer_product(mean_y, mean_y)
         st = cross_yx - outer_product(mean_y, kf%x_filt(:, t))
         model_var = rt - matmul(matmul(zt, kf%p_filt(:, :, t)), transpose(zt)) + &
            matmul(st, transpose(zt)) + matmul(zt, transpose(st))
         model_var = 0.5_dp * (model_var + transpose(model_var))
         result%var_residuals(1:n, 1:n, t) = model_var
         result%residuals(n + 1:n + m, t) = nan_value
         result%var_residuals(1:n, n + 1:n + m, t) = nan_value
         result%var_residuals(n + 1:n + m, 1:n, t) = nan_value
         result%var_residuals(n + 1:n + m, n + 1:n + m, t) = nan_value
         if (do_normalize) call normalize_observation(model, t, result%residuals(1:n, t), &
            result%var_residuals(1:n, 1:n, t))
         call standardize_residual(result%residuals(1:n, t), result%var_residuals(1:n, 1:n, t), n, missing, .false., &
            result%std_residuals(1:n, t), result%marginal_residuals(1:n, t), &
            result%block_cholesky_residuals(1:n, t), result%info)
         if (result%info /= 0) return
         result%std_residuals(n + 1:n + m, t) = nan_value
         result%marginal_residuals(n + 1:n + m, t) = nan_value
         result%block_cholesky_residuals(n + 1:n + m, t) = nan_value
         where (missing)
            result%residuals(1:n, t) = nan_value
         end where
      end do
      call finalize_full_result(model, kf, result)
   end subroutine marss_residuals_filtered

   pure subroutine marss_residuals_one_step(model, kf, result, normalize)
      type(marss_model), intent(in) :: model !! Fitted MARSS model whose residuals are conditioned on observations through time t-1.
      type(marss_kf_result), intent(in) :: kf !! Filter output supplying predictions, innovations covariance, and Kalman gains.
      type(marss_residual_result), intent(out) :: result !! Full model/state one-step-ahead residual hierarchy.
      logical, intent(in), optional :: normalize !! If true, normalize residuals using effective R and Q covariances.
      real(dp) :: at(size(model%a))
      real(dp) :: filtered_mean(size(model%z, 1))
      real(dp) :: filtered_second(size(model%z, 1), size(model%z, 1))
      real(dp) :: filtered_cross(size(model%z, 1), size(model%b, 1))
      real(dp) :: model_var(size(model%z, 1), size(model%z, 1))
      real(dp) :: predicted_mean(size(model%z, 1))
      real(dp) :: predicted_second(size(model%z, 1), size(model%z, 1))
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      real(dp) :: nan_value
      logical :: do_normalize
      logical :: missing(size(model%z, 1))
      integer :: m
      integer :: n
      integer :: t
      integer :: tt

      call initialize_full_result(model, result)
      n = size(model%y, 1)
      m = size(model%b, 1)
      tt = size(model%y, 2)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      do_normalize = .false.
      if (present(normalize)) do_normalize = normalize

      do t = 1, tt
         zt = marss_z_at(model, t)
         at = marss_a_at(model, t)
         rt = marss_r_at(model, t)
         missing = ieee_is_nan(model%y(:, t))
         call current_observation_moments(model, t, kf%x_filt(:, t), kf%p_filt(:, :, t), &
            filtered_mean, filtered_second, filtered_cross)
         predicted_mean = matmul(zt, kf%x_pred(:, t)) + at
         predicted_second = outer_product(predicted_mean, predicted_mean) + rt + &
            matmul(matmul(zt, kf%p_pred(:, :, t)), transpose(zt))
         model_var = predicted_second - outer_product(predicted_mean, predicted_mean)
         result%residuals(1:n, t) = filtered_mean - predicted_mean
         result%expected_observed_residuals(:, t) = result%residuals(1:n, t)
         result%var_observed_residuals(:, :, t) = model_var
         result%var_residuals(1:n, 1:n, t) = model_var
      end do

      do t = 1, tt
         missing = ieee_is_nan(model%y(:, t))
         if (t < tt) then
            result%residuals(n + 1:n + m, t) = matmul(kf%gain(:, :, t + 1), result%residuals(1:n, t + 1))
            result%var_residuals(n + 1:n + m, n + 1:n + m, t) = &
               matmul(matmul(kf%gain(:, :, t + 1), result%var_residuals(1:n, 1:n, t + 1)), &
               transpose(kf%gain(:, :, t + 1)))
            result%var_residuals(1:n, n + 1:n + m, t) = 0.0_dp
            result%var_residuals(n + 1:n + m, 1:n, t) = 0.0_dp
            if (do_normalize) call normalize_residual(model, t, result%residuals(:, t), result%var_residuals(:, :, t), .true.)
            call standardize_residual(result%residuals(:, t), result%var_residuals(:, :, t), n, missing, .true., &
               result%std_residuals(:, t), result%marginal_residuals(:, t), &
               result%block_cholesky_residuals(:, t), result%info)
            if (result%info /= 0) return
         else
            result%residuals(n + 1:n + m, t) = nan_value
            result%var_residuals(1:n, n + 1:n + m, t) = nan_value
            result%var_residuals(n + 1:n + m, 1:n, t) = nan_value
            result%var_residuals(n + 1:n + m, n + 1:n + m, t) = nan_value
            if (do_normalize) call normalize_observation(model, t, result%residuals(1:n, t), &
               result%var_residuals(1:n, 1:n, t))
            call standardize_residual(result%residuals(1:n, t), result%var_residuals(1:n, 1:n, t), n, missing, .false., &
               result%std_residuals(1:n, t), result%marginal_residuals(1:n, t), &
               result%block_cholesky_residuals(1:n, t), result%info)
            if (result%info /= 0) return
            result%std_residuals(:, t) = nan_value
            result%marginal_residuals(n + 1:n + m, t) = nan_value
            result%block_cholesky_residuals(n + 1:n + m, t) = nan_value
         end if
         where (missing)
            result%residuals(1:n, t) = nan_value
         end where
      end do
      call finalize_full_result(model, kf, result)
   end subroutine marss_residuals_one_step

   pure subroutine marss_residuals_harvey(model, kf, result, normalize)
      type(marss_model), intent(in) :: model !! Non-diffuse MARSS model for the Harvey disturbance-smoothing residual recursion.
      type(marss_kf_result), intent(in) :: kf !! Native MARSSkfss output supplying innovations, prediction covariances, and gains.
      type(marss_residual_result), intent(out) :: result !! Harvey tT observation/state residuals and their conditional variances.
      logical, intent(in), optional :: normalize !! If true, scale observation and state disturbances by ordered covariance factors.
      type(marss_hatyt_result) :: hat
      real(dp), allocatable :: finv(:, :)
      real(dp) :: bnext(size(model%b, 1), size(model%b, 2))
      real(dp) :: fmat(size(model%r, 1), size(model%r, 2))
      real(dp) :: jstar(size(model%b, 1), size(model%r, 1) + size(model%b, 1))
      real(dp) :: ktrans(size(model%b, 1), size(model%r, 1))
      real(dp) :: lmat(size(model%b, 1), size(model%b, 2))
      real(dp) :: nback(size(model%b, 1), size(model%b, 2), size(model%y, 2))
      real(dp) :: qnext(size(model%q, 1), size(model%q, 2))
      real(dp) :: qstar(size(model%b, 1), size(model%r, 1) + size(model%b, 1))
      real(dp) :: rback(size(model%b, 1), size(model%y, 2))
      real(dp) :: rstar(size(model%r, 1), size(model%r, 1) + size(model%b, 1))
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      real(dp) :: rwork(size(model%r, 1), size(model%r, 2))
      real(dp) :: uback(size(model%r, 1), size(model%y, 2))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      real(dp) :: zwork(size(model%z, 1), size(model%z, 2))
      real(dp) :: nan_value
      logical :: do_normalize
      logical :: missing(size(model%z, 1))
      real(dp) :: logdet
      integer :: i
      integer :: info
      integer :: m
      integer :: n
      integer :: rank
      integer :: t
      integer :: tt

      call initialize_full_result(model, result)
      if (model%diffuse) then
         result%info = 20
         return
      end if
      if (.not. kf%ok) then
         result%info = 21
         return
      end if
      n = size(model%y, 1)
      m = size(model%b, 1)
      tt = size(model%y, 2)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      do_normalize = .false.
      if (present(normalize)) do_normalize = normalize
      rback = 0.0_dp
      uback = 0.0_dp
      nback = 0.0_dp
      call marss_hatyt(model, kf, hat)

      do t = tt, 1, -1
         zt = marss_z_at(model, t)
         rt = marss_r_at(model, t)
         zwork = zt
         rwork = rt
         missing = ieee_is_nan(model%y(:, t))
         do i = 1, n
            if (missing(i)) then
               zwork(i, :) = 0.0_dp
               rwork(i, :) = 0.0_dp
               rwork(:, i) = 0.0_dp
            end if
         end do
         if (t < tt) then
            bnext = marss_b_at(model, t + 1)
            qnext = marss_q_at(model, t + 1)
         else
            bnext = marss_b_at(model, t)
            qnext = marss_q_at(model, t)
         end if
         fmat = matmul(matmul(zwork, kf%p_pred(:, :, t)), transpose(zwork)) + rwork
         fmat = 0.5_dp * (fmat + transpose(fmat))
         call psd_inverse_logdet(fmat, finv, logdet, rank, info)
         if (info /= 0) then
            result%info = 30 + t
            return
         end if
         ktrans = matmul(bnext, kf%gain(:, :, t))
         lmat = bnext - matmul(ktrans, zwork)
         rstar = 0.0_dp
         rstar(:, 1:n) = rwork
         qstar = 0.0_dp
         qstar(:, n + 1:n + m) = qnext
         jstar = qstar - matmul(ktrans, rstar)
         uback(:, t) = matmul(finv, kf%innov(:, t)) - matmul(transpose(ktrans), rback(:, t))
         if (t > 1) then
            rback(:, t - 1) = matmul(transpose(zwork), uback(:, t)) + matmul(transpose(bnext), rback(:, t))
            nback(:, :, t - 1) = matmul(matmul(transpose(zwork), finv), zwork) + &
               matmul(matmul(transpose(lmat), nback(:, :, t)), lmat)
            nback(:, :, t - 1) = 0.5_dp * (nback(:, :, t - 1) + transpose(nback(:, :, t - 1)) )
         end if
         result%residuals(:, t) = matmul(transpose(rstar), uback(:, t)) + matmul(transpose(qstar), rback(:, t))
         result%var_residuals(:, :, t) = matmul(matmul(transpose(rstar), finv), rstar) + &
            matmul(matmul(transpose(jstar), nback(:, :, t)), jstar)
         result%var_residuals(:, :, t) = 0.5_dp * &
            (result%var_residuals(:, :, t) + transpose(result%var_residuals(:, :, t)))
         if (do_normalize) then
            if (t < tt) then
               call normalize_residual(model, t, result%residuals(:, t), result%var_residuals(:, :, t), .true.)
            else
               call normalize_observation(model, t, result%residuals(1:n, t), result%var_residuals(1:n, 1:n, t))
            end if
         end if
         result%expected_observed_residuals(:, t) = result%residuals(1:n, t)
         result%var_observed_residuals(:, :, t) = hat%ot(:, :, t) - outer_product(hat%yt(:, t), hat%yt(:, t))
         call standardize_residual(result%residuals(:, t), result%var_residuals(:, :, t), n, missing, .true., &
            result%std_residuals(:, t), result%marginal_residuals(:, t), &
            result%block_cholesky_residuals(:, t), result%info)
         if (result%info /= 0) return
         do i = 1, n
            if (missing(i)) result%residuals(i, t) = nan_value
         end do
         do i = 1, n
            if (missing(i)) then
               result%var_residuals(i, 1:n, t) = nan_value
               result%var_residuals(1:n, i, t) = nan_value
            end if
         end do
      end do

      result%residuals(n + 1:n + m, tt) = nan_value
      result%var_residuals(:, n + 1:n + m, tt) = nan_value
      result%var_residuals(n + 1:n + m, :, tt) = nan_value
      result%std_residuals(:, tt) = nan_value
      result%marginal_residuals(n + 1:n + m, tt) = nan_value
      result%block_cholesky_residuals(n + 1:n + m, tt) = nan_value
      call finalize_full_result(model, kf, result)
   end subroutine marss_residuals_harvey

   pure subroutine current_observation_moments(model, t, state_mean, state_cov, mean_y, second_y, cross_yx)
      type(marss_model), intent(in) :: model !! Model defining current-time conditional observation moments.
      integer, intent(in) :: t !! One-based time at which the current observation has been assimilated.
      real(dp), intent(in) :: state_mean(:) !! State mean conditional on observations through the current time.
      real(dp), intent(in) :: state_cov(:, :) !! State covariance conditional on observations through the current time.
      real(dp), intent(out) :: mean_y(:) !! Conditional mean E[y(t)|Y(1:t)].
      real(dp), intent(out) :: second_y(:, :) !! Conditional second moment E[y(t)y(t)'|Y(1:t)].
      real(dp), intent(out) :: cross_yx(:, :) !! Conditional cross moment E[y(t)x(t)'|Y(1:t)].
      real(dp), allocatable :: correction(:, :)
      real(dp), allocatable :: r_inverse(:, :)
      real(dp), allocatable :: r_observed(:, :)
      real(dp), allocatable :: r_columns(:, :)
      real(dp) :: at(size(model%a))
      real(dp) :: covariance(size(model%z, 1), size(model%z, 1))
      real(dp) :: delta_r(size(model%z, 1), size(model%z, 1))
      real(dp) :: dz(size(model%z, 1), size(model%b, 1))
      real(dp) :: fitted_mean(size(model%z, 1))
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      real(dp) :: tolerance
      real(dp) :: y_zero(size(model%z, 1))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      logical :: missing(size(model%z, 1))
      integer, allocatable :: observed(:)
      integer :: i
      integer :: info
      integer :: j
      integer :: k
      integer :: n
      integer :: rank

      n = size(model%y, 1)
      zt = marss_z_at(model, t)
      at = marss_a_at(model, t)
      rt = marss_r_at(model, t)
      missing = ieee_is_nan(model%y(:, t))
      if (.not. any(missing)) then
         mean_y = model%y(:, t)
         second_y = outer_product(mean_y, mean_y)
         cross_yx = outer_product(mean_y, state_mean)
         return
      end if
      y_zero = 0.0_dp
      do i = 1, n
         if (.not. missing(i)) y_zero(i) = model%y(i, t)
      end do
      delta_r = 0.0_dp
      do i = 1, n
         delta_r(i, i) = 1.0_dp
      end do
      tolerance = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(rt)))
      k = 0
      do i = 1, n
         if (.not. missing(i) .and. abs(rt(i, i)) > tolerance) k = k + 1
      end do
      allocate(observed(k))
      k = 0
      do i = 1, n
         if (.not. missing(i) .and. abs(rt(i, i)) > tolerance) then
            k = k + 1
            observed(k) = i
         end if
      end do
      if (k > 0) then
         allocate(r_observed(k, k), r_columns(n, k))
         do j = 1, k
            r_columns(:, j) = rt(:, observed(j))
            do i = 1, k
               r_observed(i, j) = rt(observed(i), observed(j))
            end do
         end do
         call psd_inverse(r_observed, r_inverse, rank, info)
         if (info == 0) then
            correction = matmul(r_columns, r_inverse)
            do j = 1, k
               delta_r(:, observed(j)) = delta_r(:, observed(j)) - correction(:, j)
            end do
         end if
      end if
      fitted_mean = matmul(zt, state_mean) + at
      mean_y = y_zero - matmul(delta_r, y_zero - fitted_mean)
      dz = matmul(delta_r, zt)
      covariance = matmul(delta_r, rt) + matmul(matmul(dz, state_cov), transpose(dz))
      covariance = 0.5_dp * (covariance + transpose(covariance))
      do j = 1, n
         do i = 1, n
            if (.not. (missing(i) .and. missing(j))) covariance(i, j) = 0.0_dp
         end do
      end do
      second_y = outer_product(mean_y, mean_y) + covariance
      cross_yx = outer_product(mean_y, state_mean) + matmul(dz, state_cov)
   end subroutine current_observation_moments

   pure subroutine initialize_full_result(model, result)
      type(marss_model), intent(in) :: model !! Model supplying observation, state, and time dimensions.
      type(marss_residual_result), intent(out) :: result !! Residual result whose full arrays are allocated and initialized.
      real(dp) :: nan_value
      integer :: m
      integer :: n
      integer :: tt

      n = size(model%y, 1)
      m = size(model%b, 1)
      tt = size(model%y, 2)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      result%info = 0
      allocate(result%model_residuals(n, tt), result%state_residuals(m, tt), result%residuals(n + m, tt))
      allocate(result%var_residuals(n + m, n + m, tt), result%std_residuals(n + m, tt))
      allocate(result%marginal_residuals(n + m, tt), result%block_cholesky_residuals(n + m, tt))
      allocate(result%expected_observed_residuals(n, tt), result%var_observed_residuals(n, n, tt))
      result%model_residuals = nan_value
      result%state_residuals = nan_value
      result%residuals = nan_value
      result%var_residuals = nan_value
      result%std_residuals = nan_value
      result%marginal_residuals = nan_value
      result%block_cholesky_residuals = nan_value
      result%expected_observed_residuals = nan_value
      result%var_observed_residuals = nan_value
   end subroutine initialize_full_result

   pure subroutine finalize_full_result(model, kf, result)
      type(marss_model), intent(in) :: model !! Model whose missing observations determine legacy residual masking.
      type(marss_kf_result), intent(in) :: kf !! Filter output supplying innovation residuals and innovation covariances.
      type(marss_residual_result), intent(inout) :: result !! Full result augmented with legacy innovation/smoothed fields.
      real(dp) :: nan_value
      integer :: i
      integer :: n
      integer :: t
      integer :: tt

      n = size(model%y, 1)
      tt = size(model%y, 2)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      result%model_residuals = result%residuals(1:n, :)
      result%state_residuals = result%residuals(n + 1:, :)
      allocate(result%innovations(n, tt), result%standardized_innovations(n, tt), result%smoothed(n, tt))
      result%innovations = kf%innov
      result%standardized_innovations = nan_value
      result%smoothed = result%model_residuals
      do t = 1, tt
         do i = 1, n
            if (ieee_is_nan(model%y(i, t))) then
               result%innovations(i, t) = nan_value
            else if (kf%sigma(i, i, t) > 0.0_dp) then
               result%standardized_innovations(i, t) = kf%innov(i, t) / sqrt(kf%sigma(i, i, t))
            end if
         end do
      end do
   end subroutine finalize_full_result

   pure subroutine standardize_residual(residual, variance, nobs, missing, has_state, standard, marginal, block, info)
      real(dp), intent(in) :: residual(:) !! Residual vector before covariance standardization.
      real(dp), intent(in) :: variance(:, :) !! Residual covariance corresponding to residual.
      integer, intent(in) :: nobs !! Number of observation components at the beginning of residual.
      logical, intent(in) :: missing(:) !! Missing-observation mask for the nobs observation components.
      logical, intent(in) :: has_state !! True when state-residual components and their covariance are present.
      real(dp), intent(out) :: standard(:) !! Joint ordered-Cholesky standardized residuals, with missing entries as NaN.
      real(dp), intent(out) :: marginal(:) !! Componentwise residuals divided by marginal standard deviations.
      real(dp), intent(out) :: block(:) !! Observation joint-whitened residuals plus separately whitened state residuals.
      integer, intent(out) :: info !! Zero on success; nonzero if a residual covariance is materially indefinite.
      real(dp), allocatable :: work(:)
      real(dp), allocatable :: whitened(:)
      real(dp), allocatable :: state_whitened(:)
      real(dp) :: nan_value
      real(dp) :: scale
      real(dp) :: tol
      integer :: i
      integer :: nall
      integer :: state_info

      nall = size(residual)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(work(nall))
      work = residual
      do i = 1, min(nobs, size(missing))
         if (missing(i)) work(i) = 0.0_dp
      end do
      standard = nan_value
      marginal = nan_value
      block = nan_value
      scale = max(1.0_dp, maxval(abs(variance)))
      tol = 100.0_dp * epsilon(1.0_dp) * scale
      do i = 1, nall
         if (variance(i, i) > tol) then
            marginal(i) = work(i) / sqrt(variance(i, i))
         else if (abs(work(i)) <= sqrt(tol)) then
            marginal(i) = 0.0_dp
         end if
      end do
      call psd_cholesky_standardize(variance, work, whitened, info)
      if (info /= 0) return
      standard = whitened
      block = standard
      if (has_state .and. nall > nobs) then
         call psd_cholesky_standardize(variance(nobs + 1:nall, nobs + 1:nall), &
            work(nobs + 1:nall), state_whitened, state_info)
         if (state_info /= 0) then
            info = 10 + state_info
            return
         end if
         block(nobs + 1:nall) = state_whitened
      end if
      do i = 1, min(nobs, size(missing))
         if (missing(i)) then
            standard(i) = nan_value
            marginal(i) = nan_value
            block(i) = nan_value
         end if
      end do
      info = 0
   end subroutine standardize_residual

   pure subroutine normalize_observation(model, t, residual, variance)
      type(marss_model), intent(in) :: model !! Model supplying the effective observation covariance used for normalization.
      integer, intent(in) :: t !! One-based time whose effective R covariance is used.
      real(dp), intent(inout) :: residual(:) !! Observation residual transformed in place by R^{-1/2}.
      real(dp), intent(inout) :: variance(:, :) !! Observation residual covariance transformed in place by R^{-1/2}.
      real(dp), allocatable :: rinv_factor(:, :)
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      integer :: info
      rt = marss_r_at(model, t)
      call psd_cholesky_inverse(rt, rinv_factor, info)
      if (info /= 0) return
      residual = matmul(rinv_factor, residual)
      variance = matmul(matmul(rinv_factor, variance), transpose(rinv_factor))
   end subroutine normalize_observation

   pure subroutine normalize_residual(model, t, residual, variance, has_state)
      type(marss_model), intent(in) :: model !! Model supplying effective R(t) and Q(t+1) normalization covariances.
      integer, intent(in) :: t !! One-based residual time; the state residual uses process covariance at t+1.
      real(dp), intent(inout) :: residual(:) !! Joint observation/state residual transformed in place.
      real(dp), intent(inout) :: variance(:, :) !! Joint residual covariance transformed in place.
      logical, intent(in) :: has_state !! True when the residual includes a valid state block.
      real(dp), allocatable :: qinv_factor(:, :)
      real(dp), allocatable :: rinv_factor(:, :)
      real(dp), allocatable :: transform(:, :)
      real(dp) :: qt(size(model%q, 1), size(model%q, 2))
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      integer :: info
      integer :: m
      integer :: n
      n = size(model%y, 1)
      m = size(model%b, 1)
      allocate(transform(n + m, n + m))
      transform = 0.0_dp
      rt = marss_r_at(model, t)
      call psd_cholesky_inverse(rt, rinv_factor, info)
      if (info /= 0) return
      transform(1:n, 1:n) = rinv_factor
      if (has_state) then
         qt = marss_q_at(model, t + 1)
         call psd_cholesky_inverse(qt, qinv_factor, info)
         if (info /= 0) return
         transform(n + 1:n + m, n + 1:n + m) = qinv_factor
      end if
      residual = matmul(transform, residual)
      variance = matmul(matmul(transform, variance), transpose(transform))
   end subroutine normalize_residual

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector in the rank-one outer product.
      real(dp), intent(in) :: y(:) !! Right vector in the rank-one outer product.
      real(dp) :: a(size(x), size(y))

      a = spread(x, 2, size(y)) * spread(y, 1, size(x))
   end function outer_product

end module marss_residuals_full
