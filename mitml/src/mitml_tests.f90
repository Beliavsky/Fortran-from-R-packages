! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
! D1-D4 multiple-imputation hypothesis tests translated from mitml.
module mitml_tests
   use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
   use r_distributions, only : r_pf
   use r_kinds, only : dp
   use r_linalg, only : solve_system
   use mitml_numeric, only : covariance_columns, mean_columns, mean_cube_slices, mean_real, sample_variance, trace_matrix
   use mitml_types, only : MITML_ERR_ARGUMENT, MITML_ERR_DIMENSION, MITML_ERR_LINALG, MITML_OK, mi_test_result
   implicit none
   private

   public :: d1_test
   public :: d2_test
   public :: d3_test
   public :: d4_test
   public :: test_linear_constraints
   public :: test_transformed_constraints

contains

   subroutine d1_test(qhat, uhat, result, df_complete)
      real(dp), intent(in) :: qhat(:, :) !! Estimates of a k-dimensional estimand, shape k by m.
      real(dp), intent(in) :: uhat(:, :, :) !! Complete-data covariance estimates, shape k by k by m.
      type(mi_test_result), intent(out) :: result !! D1 F statistic, degrees of freedom, p-value, and RIV.
      real(dp), intent(in), optional :: df_complete !! Optional complete-data degrees of freedom for Reiter correction.
      real(dp), allocatable :: between(:, :)
      real(dp), allocatable :: qbar(:)
      real(dp), allocatable :: solved(:, :)
      real(dp), allocatable :: ttilde(:, :)
      real(dp), allocatable :: ubar(:, :)
      real(dp), allocatable :: x(:)
      real(dp) :: a
      real(dp) :: c0
      real(dp) :: c1
      real(dp) :: c2
      real(dp) :: r
      real(dp) :: t
      real(dp) :: vstar
      real(dp) :: z
      integer :: info
      integer :: k
      integer :: m

      call clear_test_result(result)
      k = size(qhat, 1)
      m = size(qhat, 2)
      if (k < 1 .or. m < 2) then
         call set_test_error(result, MITML_ERR_ARGUMENT, "D1 requires k >= 1 and at least two imputations")
         return
      end if
      if (size(uhat, 1) /= k .or. size(uhat, 2) /= k .or. size(uhat, 3) /= m) then
         call set_test_error(result, MITML_ERR_DIMENSION, "uhat must have shape k by k by m")
         return
      end if
      if (present(df_complete)) then
         if (df_complete <= 0.0_dp) then
            call set_test_error(result, MITML_ERR_ARGUMENT, "df_complete must be positive")
            return
         end if
      end if

      allocate(qbar(k), ubar(k, k), between(k, k), solved(k, k), ttilde(k, k), x(k))
      call mean_columns(qhat, qbar)
      call mean_cube_slices(uhat, ubar)
      call covariance_columns(qhat, between)
      call solve_system(ubar, between, solved, info)
      if (info /= 0) then
         call set_test_error(result, MITML_ERR_LINALG, "D1 could not solve with the pooled within-imputation covariance")
         return
      end if
      r = (1.0_dp + 1.0_dp / real(m, dp)) * trace_matrix(solved) / real(k, dp)
      ttilde = (1.0_dp + r) * ubar
      call solve_system(ttilde, qbar, x, info)
      if (info /= 0) then
         call set_test_error(result, MITML_ERR_LINALG, "D1 pooled covariance is singular")
         return
      end if
      result%f_value = dot_product(qbar, x) / real(k, dp)
      result%df1 = real(k, dp)
      result%riv = r
      t = real(k * (m - 1), dp)

      if (present(df_complete)) then
         if (t <= 4.0_dp) then
            call set_test_error(result, MITML_ERR_ARGUMENT, &
               "finite-df D1 correction requires k*(m-1) greater than four")
            return
         end if
         a = r * t / (t - 2.0_dp)
         vstar = ((df_complete + 1.0_dp) / (df_complete + 3.0_dp)) * df_complete
         c0 = 1.0_dp / (t - 4.0_dp)
         c1 = vstar - 2.0_dp * (1.0_dp + a)
         c2 = vstar - 4.0_dp * (1.0_dp + a)
         if (abs(c1) <= tiny(1.0_dp) .or. abs(c2) <= tiny(1.0_dp)) then
            call set_test_error(result, MITML_ERR_ARGUMENT, "finite-df D1 correction is singular for these inputs")
            return
         end if
         z = 1.0_dp / c2
         z = z + c0 * (a**2 * c1 / ((1.0_dp + a)**2 * c2))
         z = z + c0 * (8.0_dp * a**2 * c1 / ((1.0_dp + a) * c2**2) &
            + 4.0_dp * a**2 / ((1.0_dp + a) * c2))
         z = z + c0 * (4.0_dp * a**2 / (c2 * c1) + 16.0_dp * a**2 * c1 / c2**3)
         z = z + c0 * (8.0_dp * a**2 / c2**2)
         result%df2 = 4.0_dp + 1.0_dp / z
      else
         result%df2 = repeated_df(t, real(k, dp), r)
      end if
      result%p_value = r_pf(result%f_value, result%df1, result%df2, lower_tail=.false.)
      result%status = MITML_OK
      result%message = "ok"
   end subroutine d1_test

   subroutine d2_test(d, k, result, positive_riv)
      real(dp), intent(in) :: d(:) !! Complete-data Wald or likelihood-ratio chi-square statistics, one per imputation.
      integer, intent(in) :: k !! Numerator degrees of freedom, equal to the number of tested restrictions.
      type(mi_test_result), intent(out) :: result !! D2 F statistic, degrees of freedom, p-value, and RIV.
      logical, intent(in), optional :: positive_riv !! Clamp the relative increase in variance to zero when negative.
      real(dp), allocatable :: roots(:)
      real(dp) :: dbar
      real(dp) :: r
      logical :: clamp_positive
      integer :: m

      call clear_test_result(result)
      m = size(d)
      if (m < 2 .or. k < 1) then
         call set_test_error(result, MITML_ERR_ARGUMENT, "D2 requires at least two imputations and k >= 1")
         return
      end if
      if (any(d < 0.0_dp)) then
         call set_test_error(result, MITML_ERR_ARGUMENT, "D2 complete-data statistics must be nonnegative")
         return
      end if
      allocate(roots(m))
      roots = sqrt(d)
      dbar = mean_real(d)
      r = (1.0_dp + 1.0_dp / real(m, dp)) * sample_variance(roots)
      clamp_positive = .false.
      if (present(positive_riv)) clamp_positive = positive_riv
      if (clamp_positive) r = max(0.0_dp, r)
      result%f_value = (dbar / real(k, dp) - real(m + 1, dp) / real(m - 1, dp) * r) / (1.0_dp + r)
      result%df1 = real(k, dp)
      result%riv = r
      if (abs(r) <= tiny(1.0_dp)) then
         result%df2 = ieee_value(0.0_dp, ieee_positive_inf)
      else
         result%df2 = real(k, dp)**(-3.0_dp / real(m, dp)) * real(m - 1, dp) * (1.0_dp + 1.0_dp / r)**2
      end if
      result%p_value = r_pf(result%f_value, result%df1, result%df2, lower_tail=.false.)
      result%status = MITML_OK
      result%message = "ok"
   end subroutine d2_test

   subroutine d3_test(ll_model, ll_null, ll_model_pooled, ll_null_pooled, df_model, df_null, result, positive_riv)
      real(dp), intent(in) :: ll_model(:) !! Per-imputation log likelihoods for the fitted alternative model.
      real(dp), intent(in) :: ll_null(:) !! Per-imputation log likelihoods for the corresponding null model.
      real(dp), intent(in) :: ll_model_pooled(:) !! Alternative log likelihoods evaluated at pooled parameters.
      real(dp), intent(in) :: ll_null_pooled(:) !! Null log likelihoods evaluated at pooled parameters.
      integer, intent(in) :: df_model !! Number of estimated parameters in the alternative model.
      integer, intent(in) :: df_null !! Number of estimated parameters in the null model.
      type(mi_test_result), intent(out) :: result !! D3 likelihood-ratio combination result.
      logical, intent(in), optional :: positive_riv !! Clamp negative RIV values to zero when requested.
      real(dp) :: dl_bar
      real(dp) :: dl_tilde
      real(dp) :: r
      real(dp) :: t
      integer :: k
      integer :: m
      logical :: clamp_positive

      call clear_test_result(result)
      m = size(ll_model)
      k = df_model - df_null
      if (m < 2 .or. k < 1) then
         call set_test_error(result, MITML_ERR_ARGUMENT, "D3 requires at least two imputations and nested positive df difference")
         return
      end if
      if (size(ll_null) /= m .or. size(ll_model_pooled) /= m .or. size(ll_null_pooled) /= m) then
         call set_test_error(result, MITML_ERR_DIMENSION, "all D3 log-likelihood arrays must have the same length")
         return
      end if
      dl_bar = mean_real(-2.0_dp * (ll_null - ll_model))
      dl_tilde = mean_real(-2.0_dp * (ll_null_pooled - ll_model_pooled))
      r = real(m + 1, dp) / real(k * (m - 1), dp) * (dl_bar - dl_tilde)
      clamp_positive = .false.
      if (present(positive_riv)) clamp_positive = positive_riv
      if (clamp_positive) r = max(0.0_dp, r)
      result%f_value = dl_tilde / (real(k, dp) * (1.0_dp + r))
      result%df1 = real(k, dp)
      result%riv = r
      t = real(k * (m - 1), dp)
      result%df2 = repeated_df(t, real(k, dp), r)
      result%p_value = r_pf(result%f_value, result%df1, result%df2, lower_tail=.false.)
      result%status = MITML_OK
      result%message = "ok"
   end subroutine d3_test

   subroutine d4_test(ll_model, ll_null, ll_model_stacked, ll_null_stacked, df_model, df_null, result, &
      positive_riv, robust_riv)
      real(dp), intent(in) :: ll_model(:) !! Alternative-model log likelihoods fitted separately to each imputation.
      real(dp), intent(in) :: ll_null(:) !! Null-model log likelihoods fitted separately to each imputation.
      real(dp), intent(in) :: ll_model_stacked !! Alternative-model log likelihood for the stacked-data fit.
      real(dp), intent(in) :: ll_null_stacked !! Null-model log likelihood for the stacked-data fit.
      integer, intent(in) :: df_model !! Alternative-model parameter count used by D4.
      integer, intent(in) :: df_null !! Null-model parameter count used by D4.
      type(mi_test_result), intent(out) :: result !! D4 stacked-likelihood combination result.
      logical, intent(in), optional :: positive_riv !! Clamp negative nonrobust RIV values to zero.
      logical, intent(in), optional :: robust_riv !! Use Chan-Meng's robust RIV based on the alternative model likelihood.
      real(dp) :: dbar
      real(dp) :: dhat
      real(dp) :: delta_bar
      real(dp) :: delta_hat
      real(dp) :: r
      integer :: h
      integer :: k
      integer :: m
      logical :: clamp_positive
      logical :: robust

      call clear_test_result(result)
      m = size(ll_model)
      k = df_model - df_null
      h = df_model
      if (m < 2 .or. k < 1 .or. h < 1) then
         call set_test_error(result, MITML_ERR_ARGUMENT, "D4 requires at least two imputations and a nested model")
         return
      end if
      if (size(ll_null) /= m) then
         call set_test_error(result, MITML_ERR_DIMENSION, "D4 log-likelihood arrays must have equal length")
         return
      end if
      dbar = mean_real(-2.0_dp * (ll_null - ll_model))
      dhat = -2.0_dp * (ll_null_stacked - ll_model_stacked)
      robust = .false.
      if (present(robust_riv)) robust = robust_riv
      clamp_positive = .false.
      if (present(positive_riv)) clamp_positive = positive_riv

      if (robust) then
         delta_bar = 2.0_dp * mean_real(ll_model)
         delta_hat = 2.0_dp * ll_model_stacked
         r = real(m + 1, dp) / real(h * (m - 1), dp) * (delta_bar - delta_hat)
         result%df2 = repeated_simple_df(real(h * (m - 1), dp), r)
      else
         r = real(m + 1, dp) / real(k * (m - 1), dp) * (dbar - dhat)
         if (clamp_positive) r = max(0.0_dp, r)
         result%df2 = repeated_simple_df(real(k * (m - 1), dp), r)
      end if
      result%f_value = dhat / (real(k, dp) * (1.0_dp + r))
      result%df1 = real(k, dp)
      result%riv = r
      result%p_value = r_pf(result%f_value, result%df1, result%df2, lower_tail=.false.)
      result%status = MITML_OK
      result%message = "ok"
   end subroutine d4_test

   subroutine test_linear_constraints(qhat, uhat, a, b, method, result, df_complete, positive_riv)
      real(dp), intent(in) :: qhat(:, :) !! Original parameter estimates, shape p by m.
      real(dp), intent(in) :: uhat(:, :, :) !! Original covariance matrices, shape p by p by m.
      real(dp), intent(in) :: a(:, :) !! Linear restriction matrix, shape k by p.
      real(dp), intent(in) :: b(:) !! Restriction target vector; tests A*theta - b equal to zero.
      integer, intent(in) :: method !! Pooling method: 1 for D1 or 2 for D2.
      type(mi_test_result), intent(out) :: result !! Multiple-imputation test of the linear restrictions.
      real(dp), intent(in), optional :: df_complete !! Complete-data degrees of freedom used only by D1.
      logical, intent(in), optional :: positive_riv !! D2 option to clamp a negative RIV to zero.
      real(dp), allocatable :: transformed_q(:, :)
      real(dp), allocatable :: transformed_u(:, :, :)
      integer :: i
      integer :: k
      integer :: m
      integer :: p

      call clear_test_result(result)
      p = size(qhat, 1)
      m = size(qhat, 2)
      k = size(a, 1)
      if (size(a, 2) /= p .or. size(b) /= k) then
         call set_test_error(result, MITML_ERR_DIMENSION, "constraint matrix/vector shapes are incompatible with qhat")
         return
      end if
      if (size(uhat, 1) /= p .or. size(uhat, 2) /= p .or. size(uhat, 3) /= m) then
         call set_test_error(result, MITML_ERR_DIMENSION, "uhat shape is incompatible with qhat")
         return
      end if
      allocate(transformed_q(k, m), transformed_u(k, k, m))
      do i = 1, m
         transformed_q(:, i) = matmul(a, qhat(:, i)) - b
         transformed_u(:, :, i) = matmul(a, matmul(uhat(:, :, i), transpose(a)))
      end do
      call test_transformed_constraints(transformed_q, transformed_u, method, result, df_complete, positive_riv)
   end subroutine test_linear_constraints

   subroutine test_transformed_constraints(qhat, uhat, method, result, df_complete, positive_riv)
      real(dp), intent(in) :: qhat(:, :) !! Already transformed constraint estimates, shape k by m.
      real(dp), intent(in) :: uhat(:, :, :) !! Delta-method covariance matrices, shape k by k by m.
      integer, intent(in) :: method !! Pooling method: 1 for D1 or 2 for D2.
      type(mi_test_result), intent(out) :: result !! Test result for the transformed constraints.
      real(dp), intent(in), optional :: df_complete !! Complete-data degrees of freedom for D1.
      logical, intent(in), optional :: positive_riv !! D2 option to clamp the RIV at zero.
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: x(:)
      integer :: i
      integer :: info
      integer :: k
      integer :: m

      k = size(qhat, 1)
      m = size(qhat, 2)
      if (method == 1) then
         if (present(df_complete)) then
            call d1_test(qhat, uhat, result, df_complete)
         else
            call d1_test(qhat, uhat, result)
         end if
         return
      end if
      if (method /= 2) then
         call clear_test_result(result)
         call set_test_error(result, MITML_ERR_ARGUMENT, "constraint method must be 1 (D1) or 2 (D2)")
         return
      end if
      if (size(uhat, 1) /= k .or. size(uhat, 2) /= k .or. size(uhat, 3) /= m) then
         call clear_test_result(result)
         call set_test_error(result, MITML_ERR_DIMENSION, "transformed uhat must have shape k by k by m")
         return
      end if
      allocate(d(m), x(k))
      do i = 1, m
         call solve_system(uhat(:, :, i), qhat(:, i), x, info)
         if (info /= 0) then
            call clear_test_result(result)
            call set_test_error(result, MITML_ERR_LINALG, "D2 constraint covariance is singular")
            return
         end if
         d(i) = dot_product(qhat(:, i), x)
      end do
      if (present(positive_riv)) then
         call d2_test(d, k, result, positive_riv)
      else
         call d2_test(d, k, result)
      end if
   end subroutine test_transformed_constraints

   pure real(dp) function repeated_df(t, k, r) result(value)
      real(dp), intent(in) :: t !! Product k*(m-1) in the Li-Meng repeated-imputation df formula.
      real(dp), intent(in) :: k !! Numerator dimension used by the small-t branch.
      real(dp), intent(in) :: r !! Relative increase in variance.

      if (abs(r) <= tiny(1.0_dp)) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
      else if (t > 4.0_dp) then
         value = 4.0_dp + (t - 4.0_dp) * (1.0_dp + (1.0_dp - 2.0_dp / t) / r)**2
      else
         value = t * (1.0_dp + 1.0_dp / k) * (1.0_dp + 1.0_dp / r)**2 / 2.0_dp
      end if
   end function repeated_df

   pure real(dp) function repeated_simple_df(t, r) result(value)
      real(dp), intent(in) :: t !! Product of model dimension and m - 1 in the D4 df formula.
      real(dp), intent(in) :: r !! Relative increase in variance.

      if (abs(r) <= tiny(1.0_dp)) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
      else
         value = t * (1.0_dp + 1.0_dp / r)**2
      end if
   end function repeated_simple_df

   subroutine clear_test_result(result)
      type(mi_test_result), intent(out) :: result !! Test result reset before a calculation.

      result%f_value = 0.0_dp
      result%df1 = 0.0_dp
      result%df2 = 0.0_dp
      result%p_value = 1.0_dp
      result%riv = 0.0_dp
      result%status = MITML_OK
      result%message = ""
   end subroutine clear_test_result

   subroutine set_test_error(result, status, message)
      type(mi_test_result), intent(inout) :: result !! Test result receiving an error status and message.
      integer, intent(in) :: status !! MITML status code describing the failure.
      character(len=*), intent(in) :: message !! Human-readable explanation of the failure.

      result%status = status
      result%message = message
      result%f_value = ieee_value(0.0_dp, ieee_quiet_nan)
      result%df1 = ieee_value(0.0_dp, ieee_quiet_nan)
      result%df2 = ieee_value(0.0_dp, ieee_quiet_nan)
      result%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
      result%riv = ieee_value(0.0_dp, ieee_quiet_nan)
   end subroutine set_test_error

end module mitml_tests
