! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package roll by Jason Foster; see PROVENANCE.md and NOTICE.md.
module roll_multivariate
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use r_linalg, only : solve_system
   use roll_univariate, only : roll_pair_vec
   implicit none
   private

   type, public :: roll_lm_result
      real(dp), allocatable :: coefficients(:, :, :) !! `(row, coefficient, response)` rolling estimates.
      real(dp), allocatable :: r_squared(:, :)       !! `(row, response)` rolling coefficients of determination.
      real(dp), allocatable :: std_error(:, :, :)    !! `(row, coefficient, response)` rolling standard errors.
   end type roll_lm_result

   public :: roll_cov_vec, roll_cor_vec, roll_crossprod_vec
   public :: roll_cov_mat, roll_cor_mat, roll_crossprod_mat
   public :: roll_lm_vec, roll_lm_mat

contains

   pure function na_real() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function na_real

   pure subroutine roll_cov_vec(x, width, result, y, weights, center, scale, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! First numeric time series.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Rolling covariance or correlation when `scale` is true.
      real(dp), intent(in), optional :: y(:) !! Second series; when absent `x` is paired with itself.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: center !! Subtract weighted means when true; defaults to true.
      logical, intent(in), optional :: scale !! Normalize to correlation when true; defaults to false.
      integer, intent(in), optional :: min_obs !! Minimum pairwise-complete observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vector pairs are intrinsically pairwise.
      logical, intent(in), optional :: na_restore !! Restore a missing current pair when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      real(dp), allocatable :: yy(:)
      logical :: use_center, use_scale

      use_center = .true.
      if (present(center)) use_center = center
      use_scale = .false.
      if (present(scale)) use_scale = scale
      if (present(y)) then
         allocate(yy(size(y)))
         yy = y
      else
         allocate(yy(size(x)))
         yy = x
      end if
      call roll_pair_vec(x, yy, width, result, weights, use_center, use_scale, min_obs, na_restore, .false.)
   end subroutine roll_cov_vec

   pure subroutine roll_cor_vec(x, width, result, y, weights, center, scale, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! First numeric time series.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Rolling correlations.
      real(dp), intent(in), optional :: y(:) !! Second series; when absent `x` is paired with itself.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: center !! Subtract weighted means when true; defaults to true.
      logical, intent(in), optional :: scale !! Accepted for R API parity; correlation always scales.
      integer, intent(in), optional :: min_obs !! Minimum pairwise-complete observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vector pairs are intrinsically pairwise.
      logical, intent(in), optional :: na_restore !! Restore a missing current pair when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      logical :: use_center

      use_center = .true.
      if (present(center)) use_center = center
      call roll_cov_vec(x, width, result, y, weights, use_center, .true., min_obs, complete_obs, na_restore, online)
   end subroutine roll_cor_vec

   pure subroutine roll_crossprod_vec(x, width, result, y, weights, center, scale, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! First numeric time series.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Rolling centered/scaled crossproducts.
      real(dp), intent(in), optional :: y(:) !! Second series; when absent `x` is paired with itself.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: center !! Subtract weighted means when true; defaults to false.
      logical, intent(in), optional :: scale !! Normalize to a correlation-like crossproduct when true.
      integer, intent(in), optional :: min_obs !! Minimum pairwise-complete observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vector pairs are intrinsically pairwise.
      logical, intent(in), optional :: na_restore !! Restore a missing current pair when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      real(dp), allocatable :: yy(:)
      logical :: use_center, use_scale

      use_center = .false.
      if (present(center)) use_center = center
      use_scale = .false.
      if (present(scale)) use_scale = scale
      if (present(y)) then
         allocate(yy(size(y)))
         yy = y
      else
         allocate(yy(size(x)))
         yy = x
      end if
      call roll_pair_vec(x, yy, width, result, weights, use_center, use_scale, min_obs, na_restore, .true.)
   end subroutine roll_crossprod_vec

   pure subroutine roll_cov_mat(x, width, result, y, weights, center, scale, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! First numeric matrix with observations in rows.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :, :) !! `(x-column, y-column, row)` rolling covariance cube.
      real(dp), intent(in), optional :: y(:, :) !! Second matrix; when absent `x` is paired with itself.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: center !! Subtract weighted means when true; defaults to true.
      logical, intent(in), optional :: scale !! Normalize to correlations when true; defaults to false.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude a row for all pairs when any source value is missing.
      logical, intent(in), optional :: na_restore !! Restore missing current pair cells when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      logical :: use_center, use_scale

      use_center = .true.
      if (present(center)) use_center = center
      use_scale = .false.
      if (present(scale)) use_scale = scale
      call pair_matrix_core(x, width, result, y, weights, use_center, use_scale, min_obs, complete_obs, na_restore, .false.)
   end subroutine roll_cov_mat

   pure subroutine roll_cor_mat(x, width, result, y, weights, center, scale, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! First numeric matrix with observations in rows.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :, :) !! `(x-column, y-column, row)` rolling correlation cube.
      real(dp), intent(in), optional :: y(:, :) !! Second matrix; when absent `x` is paired with itself.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: center !! Subtract weighted means when true; defaults to true.
      logical, intent(in), optional :: scale !! Accepted for R API parity; correlation always scales.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude a row for all pairs when any source value is missing.
      logical, intent(in), optional :: na_restore !! Restore missing current pair cells when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      logical :: use_center

      use_center = .true.
      if (present(center)) use_center = center
      call pair_matrix_core(x, width, result, y, weights, use_center, .true., min_obs, complete_obs, na_restore, .false.)
   end subroutine roll_cor_mat

   pure subroutine roll_crossprod_mat(x, width, result, y, weights, center, scale, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! First numeric matrix with observations in rows.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :, :) !! `(x-column, y-column, row)` rolling crossproduct cube.
      real(dp), intent(in), optional :: y(:, :) !! Second matrix; when absent `x` is paired with itself.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: center !! Subtract weighted means when true; defaults to false.
      logical, intent(in), optional :: scale !! Normalize to correlation-like values when true.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude a row for all pairs when any source value is missing.
      logical, intent(in), optional :: na_restore !! Restore missing current pair cells when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      logical :: use_center, use_scale

      use_center = .false.
      if (present(center)) use_center = center
      use_scale = .false.
      if (present(scale)) use_scale = scale
      call pair_matrix_core(x, width, result, y, weights, use_center, use_scale, min_obs, complete_obs, na_restore, .true.)
   end subroutine roll_crossprod_mat

   pure subroutine pair_matrix_core(x, width, result, y, weights, center, scale, min_obs, complete_obs, na_restore, cross_mode)
      real(dp), intent(in) :: x(:, :) !! First numeric matrix with observations in rows.
      integer, intent(in) :: width !! Rolling window width.
      real(dp), allocatable, intent(out) :: result(:, :, :) !! Pairwise rolling cube.
      real(dp), intent(in), optional :: y(:, :) !! Optional second matrix; absent means `x`.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in) :: center !! Subtract weighted means when true.
      logical, intent(in) :: scale !! Normalize pairwise products when true.
      integer, intent(in), optional :: min_obs !! Minimum usable observations.
      logical, intent(in), optional :: complete_obs !! Exclude rows with any source missing when true.
      logical, intent(in), optional :: na_restore !! Restore current-pair NaNs when true.
      logical, intent(in) :: cross_mode !! Return crossproducts instead of unbiased covariances when true.
      real(dp), allocatable :: yy(:, :), xv(:), yv(:), temp(:)
      integer :: i, j, k, ny
      logical :: complete, restore, bad_row

      if (present(y)) then
         allocate(yy(size(y, 1), size(y, 2)))
         yy = y
      else
         allocate(yy(size(x, 1), size(x, 2)))
         yy = x
      end if
      ny = size(yy, 2)
      complete = .true.
      if (present(complete_obs)) complete = complete_obs
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x, 2), ny, size(x, 1)), xv(size(x, 1)), yv(size(x, 1)))
      do j = 1, size(x, 2)
         do k = 1, ny
            xv = x(:, j)
            yv = yy(:, k)
            if (complete) then
               do i = 1, size(x, 1)
                  bad_row = any(ieee_is_nan(x(i, :))) .or. any(ieee_is_nan(yy(i, :)))
                  if (bad_row) then
                     xv(i) = na_real()
                     yv(i) = na_real()
                  end if
               end do
            end if
            call roll_pair_vec(xv, yv, width, temp, weights, center, scale, min_obs, .false., cross_mode)
            result(j, k, :) = temp
            if (restore) then
               do i = 1, size(x, 1)
                  if (ieee_is_nan(x(i, j)) .or. ieee_is_nan(yy(i, k))) result(j, k, i) = na_real()
               end do
            end if
         end do
      end do
   end subroutine pair_matrix_core

   pure subroutine roll_lm_vec(x, y, width, fit, weights, intercept, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Single predictor series.
      real(dp), intent(in) :: y(:) !! Single response series with the same length as `x`.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      type(roll_lm_result), intent(out) :: fit !! Rolling coefficients, R-squared values, and standard errors.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: intercept !! Include an intercept when true; defaults to true.
      integer, intent(in), optional :: min_obs !! Minimum complete observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for parity; upstream always uses complete rows for `roll_lm`.
      logical, intent(in), optional :: na_restore !! Restore current-row missing values when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      real(dp), allocatable :: xx(:, :), yy(:, :)

      allocate(xx(size(x), 1), yy(size(y), 1))
      xx(:, 1) = x
      yy(:, 1) = y
      call roll_lm_mat(xx, yy, width, fit, weights, intercept, min_obs, complete_obs, na_restore, online)
   end subroutine roll_lm_vec

   pure subroutine roll_lm_mat(x, y, width, fit, weights, intercept, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with observations in rows.
      real(dp), intent(in) :: y(:, :) !! Response matrix with the same row count as `x`.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      type(roll_lm_result), intent(out) :: fit !! Rolling coefficients, R-squared values, and standard errors.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: intercept !! Include an intercept when true; defaults to true.
      integer, intent(in), optional :: min_obs !! Minimum complete observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for parity; upstream always uses complete rows for `roll_lm`.
      logical, intent(in), optional :: na_restore !! Restore current-row missing values when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      integer :: i, j, k, c, p, q, n_obs, min_n, df_resid, info
      real(dp) :: sum_w, w, mean_y, yy, rsq, var_resid, quad
      real(dp), allocatable :: mean_x(:), a(:, :), b(:), beta(:), ainv(:, :), ident(:, :)
      logical :: use_intercept, restore, valid, current_missing

      use_intercept = .true.
      if (present(intercept)) use_intercept = intercept
      restore = .false.
      if (present(na_restore)) restore = na_restore
      min_n = width
      if (present(min_obs)) min_n = min_obs
      p = size(x, 2)
      q = p + merge(1, 0, use_intercept)
      allocate(fit%coefficients(size(x, 1), q, size(y, 2)))
      allocate(fit%r_squared(size(x, 1), size(y, 2)))
      allocate(fit%std_error(size(x, 1), q, size(y, 2)))
      fit%coefficients = na_real()
      fit%r_squared = na_real()
      fit%std_error = na_real()
      allocate(mean_x(p), a(p, p), b(p), beta(p), ainv(p, p), ident(p, p))
      ident = 0.0_dp
      do j = 1, p
         ident(j, j) = 1.0_dp
      end do

      do k = 1, size(y, 2)
         do i = 1, size(x, 1)
            current_missing = ieee_is_nan(y(i, k)) .or. any(ieee_is_nan(x(i, :)))
            if (restore .and. current_missing) cycle
            sum_w = 0.0_dp
            mean_x = 0.0_dp
            mean_y = 0.0_dp
            n_obs = 0
            do c = 0, min(width - 1, i - 1)
               valid = .not. ieee_is_nan(y(i - c, k)) .and. .not. any(ieee_is_nan(x(i - c, :)))
               if (valid) then
                  w = 1.0_dp
                  if (present(weights)) w = weights(size(weights) - c)
                  sum_w = sum_w + w
                  mean_x = mean_x + w * x(i - c, :)
                  mean_y = mean_y + w * y(i - c, k)
                  n_obs = n_obs + 1
               end if
            end do
            if (n_obs < min_n .or. n_obs < q .or. sum_w == 0.0_dp) cycle
            if (use_intercept) then
               mean_x = mean_x / sum_w
               mean_y = mean_y / sum_w
            else
               mean_x = 0.0_dp
               mean_y = 0.0_dp
            end if
            a = 0.0_dp
            b = 0.0_dp
            yy = 0.0_dp
            do c = 0, min(width - 1, i - 1)
               valid = .not. ieee_is_nan(y(i - c, k)) .and. .not. any(ieee_is_nan(x(i - c, :)))
               if (valid) then
                  w = 1.0_dp
                  if (present(weights)) w = weights(size(weights) - c)
                  do j = 1, p
                     b(j) = b(j) + w * (x(i - c, j) - mean_x(j)) * (y(i - c, k) - mean_y)
                     a(j, :) = a(j, :) + w * (x(i - c, j) - mean_x(j)) * (x(i - c, :) - mean_x)
                  end do
                  yy = yy + w * (y(i - c, k) - mean_y) ** 2
               end if
            end do
            call solve_system(a, b, beta, info)
            if (info /= 0) cycle
            if (use_intercept) then
               fit%coefficients(i, 1, k) = mean_y - dot_product(mean_x, beta)
               fit%coefficients(i, 2:q, k) = beta
            else
               fit%coefficients(i, 1:q, k) = beta
            end if
            if (yy > epsilon(1.0_dp)) then
               quad = dot_product(beta, matmul(a, beta))
               rsq = quad / yy
               fit%r_squared(i, k) = rsq
            else
               cycle
            end if
            df_resid = n_obs - q
            if (df_resid <= 0) cycle
            call solve_system(a, ident, ainv, info)
            if (info /= 0) cycle
            var_resid = (1.0_dp - rsq) * yy / real(df_resid, dp)
            if (var_resid < 0.0_dp .and. abs(var_resid) <= epsilon(1.0_dp)) var_resid = 0.0_dp
            if (var_resid < 0.0_dp) cycle
            if (use_intercept) then
               fit%std_error(i, 1, k) = sqrt(var_resid * (1.0_dp / sum_w + dot_product(mean_x, matmul(ainv, mean_x))))
               do j = 1, p
                  fit%std_error(i, j + 1, k) = sqrt(max(0.0_dp, var_resid * ainv(j, j)))
               end do
            else
               do j = 1, p
                  fit%std_error(i, j, k) = sqrt(max(0.0_dp, var_resid * ainv(j, j)))
               end do
            end if
         end do
      end do
   end subroutine roll_lm_mat

end module roll_multivariate
