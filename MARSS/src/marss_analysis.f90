! SPDX-License-Identifier: GPL-2.0-only
module marss_analysis
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_kf_result, marss_residual_result, marss_hatyt_result, marss_constraints
   use marss_utils, only : count_observed, normal_quantile
   use marss_covariance, only : psd_inverse
   use marss_kalman, only : marss_kfss
   use marss_parameters, only : marss_z_at, marss_a_at, marss_r_at, marss_v0_effective
   use r_linalg, only : inverse_matrix
   use marss_constraints_mod, only : marss_constraints_parameter_count, marss_vectorize_free, marss_unvectorize_free
   implicit none
   private
   public :: marss_aic
   public :: marss_hatyt
   public :: marss_residuals
   public :: marss_vectorizeparam
   public :: marss_unvectorizeparam
   public :: marss_hessian
   public :: marss_fisher_i
   public :: marss_param_cis
   public :: marss_parameter_count

contains

   pure elemental subroutine marss_aic(loglik, nparam, nsample, aic, aicc)
      real(dp), intent(in) :: loglik !! Maximized log likelihood.
      integer, intent(in) :: nparam !! Number of estimated parameters K.
      integer, intent(in) :: nsample !! Number of nonmissing scalar observations.
      real(dp), intent(out) :: aic !! Akaike information criterion, -2 logLik + 2 K.
      real(dp), intent(out) :: aicc !! Small-sample AIC correction, huge if nsample <= K+1.

      aic = -2.0_dp * loglik + 2.0_dp * real(nparam, dp)
      if (nsample > nparam + 1) then
         aicc = -2.0_dp * loglik + 2.0_dp * real(nparam, dp) * &
            real(nsample, dp) / real(nsample - nparam - 1, dp)
      else
         aicc = huge(1.0_dp)
      end if
   end subroutine marss_aic

   pure integer function marss_parameter_count(model, constraints) result(nparam)
      type(marss_model), intent(in) :: model !! Model whose active numerical or affine-free parameterization is counted.
      type(marss_constraints), intent(in), optional :: constraints !! Optional affine structure for beta-coordinate counting.
      integer :: m
      integer :: n
      integer :: tt

      if (present(constraints)) then
         nparam = marss_constraints_parameter_count(constraints)
         return
      end if
      m = size(model%b, 1)
      n = size(model%z, 1)
      tt = size(model%y, 2)
      nparam = 0
      if (allocated(model%b_t)) then
         nparam = nparam + size(model%b_t)
      else
         nparam = nparam + size(model%b)
      end if
      if (allocated(model%u_t)) then
         nparam = nparam + size(model%u_t)
      else
         nparam = nparam + size(model%u)
      end if
      if (allocated(model%q_noise)) then
         if (allocated(model%q_noise_t)) then
            nparam = nparam + tt * size(model%q_noise, 1) * (size(model%q_noise, 1) + 1) / 2
         else
            nparam = nparam + size(model%q_noise, 1) * (size(model%q_noise, 1) + 1) / 2
         end if
      else if (allocated(model%q_t)) then
         nparam = nparam + tt * m * (m + 1) / 2
      else
         nparam = nparam + m * (m + 1) / 2
      end if
      if (allocated(model%z_t)) then
         nparam = nparam + size(model%z_t)
      else
         nparam = nparam + size(model%z)
      end if
      if (allocated(model%a_t)) then
         nparam = nparam + size(model%a_t)
      else
         nparam = nparam + size(model%a)
      end if
      if (allocated(model%r_noise)) then
         if (allocated(model%r_noise_t)) then
            nparam = nparam + tt * size(model%r_noise, 1) * (size(model%r_noise, 1) + 1) / 2
         else
            nparam = nparam + size(model%r_noise, 1) * (size(model%r_noise, 1) + 1) / 2
         end if
      else if (allocated(model%r_t)) then
         nparam = nparam + tt * n * (n + 1) / 2
      else
         nparam = nparam + n * (n + 1) / 2
      end if
      nparam = nparam + size(model%x0)
      if (allocated(model%v0_noise)) then
         nparam = nparam + size(model%v0_noise, 1) * (size(model%v0_noise, 1) + 1) / 2
      else
         nparam = nparam + m * (m + 1) / 2
      end if
      if (allocated(model%c_coef_t)) then
         nparam = nparam + size(model%c_coef_t)
      else if (allocated(model%c_coef)) then
         nparam = nparam + size(model%c_coef)
      end if
      if (allocated(model%d_coef_t)) then
         nparam = nparam + size(model%d_coef_t)
      else if (allocated(model%d_coef)) then
         nparam = nparam + size(model%d_coef)
      end if
      if (allocated(model%g_t)) then
         nparam = nparam + size(model%g_t)
      else if (allocated(model%g)) then
         nparam = nparam + size(model%g)
      end if
      if (allocated(model%h_t)) then
         nparam = nparam + size(model%h_t)
      else if (allocated(model%h)) then
         nparam = nparam + size(model%h)
      end if
      if (allocated(model%l)) nparam = nparam + size(model%l)
   end function marss_parameter_count

   pure subroutine marss_vectorizeparam(model, parvec, constraints, info)
      type(marss_model), intent(in) :: model !! Model whose active numerical blocks or affine free coordinates are packed.
      real(dp), allocatable, intent(out) :: parvec(:) !! Packed raw values, or beta coordinates when constraints are supplied.
      type(marss_constraints), intent(in), optional :: constraints !! Optional f+D*beta structure selecting free coordinates only.
      integer, intent(out), optional :: info !! Zero on success; nonzero if the model violates the supplied affine constraints.
      integer :: local_info
      integer :: k
      integer :: t

      if (present(constraints)) then
         call marss_vectorize_free(model, constraints, parvec, local_info)
         if (present(info)) info = local_info
         return
      end if
      allocate(parvec(marss_parameter_count(model)))
      k = 0
      if (allocated(model%b_t)) then
         do t = 1, size(model%b_t, 3)
            call pack_full(model%b_t(:, :, t), parvec, k)
         end do
      else
         call pack_full(model%b, parvec, k)
      end if
      if (allocated(model%u_t)) then
         do t = 1, size(model%u_t, 2)
            call pack_vector(model%u_t(:, t), parvec, k)
         end do
      else
         call pack_vector(model%u, parvec, k)
      end if
      if (allocated(model%q_noise)) then
         if (allocated(model%q_noise_t)) then
            do t = 1, size(model%q_noise_t, 3)
               call pack_lower(model%q_noise_t(:, :, t), parvec, k)
            end do
         else
            call pack_lower(model%q_noise, parvec, k)
         end if
      else if (allocated(model%q_t)) then
         do t = 1, size(model%q_t, 3)
            call pack_lower(model%q_t(:, :, t), parvec, k)
         end do
      else
         call pack_lower(model%q, parvec, k)
      end if
      if (allocated(model%z_t)) then
         do t = 1, size(model%z_t, 3)
            call pack_full(model%z_t(:, :, t), parvec, k)
         end do
      else
         call pack_full(model%z, parvec, k)
      end if
      if (allocated(model%a_t)) then
         do t = 1, size(model%a_t, 2)
            call pack_vector(model%a_t(:, t), parvec, k)
         end do
      else
         call pack_vector(model%a, parvec, k)
      end if
      if (allocated(model%r_noise)) then
         if (allocated(model%r_noise_t)) then
            do t = 1, size(model%r_noise_t, 3)
               call pack_lower(model%r_noise_t(:, :, t), parvec, k)
            end do
         else
            call pack_lower(model%r_noise, parvec, k)
         end if
      else if (allocated(model%r_t)) then
         do t = 1, size(model%r_t, 3)
            call pack_lower(model%r_t(:, :, t), parvec, k)
         end do
      else
         call pack_lower(model%r, parvec, k)
      end if
      call pack_vector(model%x0, parvec, k)
      if (allocated(model%v0_noise)) then
         call pack_lower(model%v0_noise, parvec, k)
      else
         call pack_lower(model%v0, parvec, k)
      end if
      if (allocated(model%c_coef_t)) then
         do t = 1, size(model%c_coef_t, 3)
            call pack_full(model%c_coef_t(:, :, t), parvec, k)
         end do
      else if (allocated(model%c_coef)) then
         call pack_full(model%c_coef, parvec, k)
      end if
      if (allocated(model%d_coef_t)) then
         do t = 1, size(model%d_coef_t, 3)
            call pack_full(model%d_coef_t(:, :, t), parvec, k)
         end do
      else if (allocated(model%d_coef)) then
         call pack_full(model%d_coef, parvec, k)
      end if
      if (allocated(model%g_t)) then
         do t = 1, size(model%g_t, 3)
            call pack_full(model%g_t(:, :, t), parvec, k)
         end do
      else if (allocated(model%g)) then
         call pack_full(model%g, parvec, k)
      end if
      if (allocated(model%h_t)) then
         do t = 1, size(model%h_t, 3)
            call pack_full(model%h_t(:, :, t), parvec, k)
         end do
      else if (allocated(model%h)) then
         call pack_full(model%h, parvec, k)
      end if
      if (allocated(model%l)) call pack_full(model%l, parvec, k)
      if (present(info)) info = 0
   end subroutine marss_vectorizeparam

   pure subroutine marss_unvectorizeparam(template, parvec, model, info, constraints)
      type(marss_model), intent(in) :: template !! Template supplying data, active layout, dimensions, and fixed structure.
      real(dp), intent(in) :: parvec(:) !! Raw packed vector, or beta coordinates when constraints are supplied.
      type(marss_model), intent(out) :: model !! Model reconstructed from the supplied coordinate vector.
      integer, intent(out) :: info !! Zero on success, nonzero when vector length or resulting model structure is inconsistent.
      type(marss_constraints), intent(in), optional :: constraints !! Optional affine structure for beta-coordinate unpacking.
      integer :: k
      integer :: t

      if (present(constraints)) then
         if (size(parvec) /= marss_constraints_parameter_count(constraints)) then
            info = 1
            return
         end if
         call marss_unvectorize_free(template, constraints, parvec, model, info)
         return
      end if
      model = template
      if (size(parvec) /= marss_parameter_count(template)) then
         info = 1
         return
      end if
      k = 0
      if (allocated(model%b_t)) then
         do t = 1, size(model%b_t, 3)
            call unpack_full(parvec, k, model%b_t(:, :, t))
         end do
      else
         call unpack_full(parvec, k, model%b)
      end if
      if (allocated(model%u_t)) then
         do t = 1, size(model%u_t, 2)
            call unpack_vector(parvec, k, model%u_t(:, t))
         end do
      else
         call unpack_vector(parvec, k, model%u)
      end if
      if (allocated(model%q_noise)) then
         if (allocated(model%q_noise_t)) then
            do t = 1, size(model%q_noise_t, 3)
               call unpack_lower(parvec, k, model%q_noise_t(:, :, t))
            end do
         else
            call unpack_lower(parvec, k, model%q_noise)
         end if
      else if (allocated(model%q_t)) then
         do t = 1, size(model%q_t, 3)
            call unpack_lower(parvec, k, model%q_t(:, :, t))
         end do
      else
         call unpack_lower(parvec, k, model%q)
      end if
      if (allocated(model%z_t)) then
         do t = 1, size(model%z_t, 3)
            call unpack_full(parvec, k, model%z_t(:, :, t))
         end do
      else
         call unpack_full(parvec, k, model%z)
      end if
      if (allocated(model%a_t)) then
         do t = 1, size(model%a_t, 2)
            call unpack_vector(parvec, k, model%a_t(:, t))
         end do
      else
         call unpack_vector(parvec, k, model%a)
      end if
      if (allocated(model%r_noise)) then
         if (allocated(model%r_noise_t)) then
            do t = 1, size(model%r_noise_t, 3)
               call unpack_lower(parvec, k, model%r_noise_t(:, :, t))
            end do
         else
            call unpack_lower(parvec, k, model%r_noise)
         end if
      else if (allocated(model%r_t)) then
         do t = 1, size(model%r_t, 3)
            call unpack_lower(parvec, k, model%r_t(:, :, t))
         end do
      else
         call unpack_lower(parvec, k, model%r)
      end if
      call unpack_vector(parvec, k, model%x0)
      if (allocated(model%v0_noise)) then
         call unpack_lower(parvec, k, model%v0_noise)
      else
         call unpack_lower(parvec, k, model%v0)
      end if
      if (allocated(model%c_coef_t)) then
         do t = 1, size(model%c_coef_t, 3)
            call unpack_full(parvec, k, model%c_coef_t(:, :, t))
         end do
      else if (allocated(model%c_coef)) then
         call unpack_full(parvec, k, model%c_coef)
      end if
      if (allocated(model%d_coef_t)) then
         do t = 1, size(model%d_coef_t, 3)
            call unpack_full(parvec, k, model%d_coef_t(:, :, t))
         end do
      else if (allocated(model%d_coef)) then
         call unpack_full(parvec, k, model%d_coef)
      end if
      if (allocated(model%g_t)) then
         do t = 1, size(model%g_t, 3)
            call unpack_full(parvec, k, model%g_t(:, :, t))
         end do
      else if (allocated(model%g)) then
         call unpack_full(parvec, k, model%g)
      end if
      if (allocated(model%h_t)) then
         do t = 1, size(model%h_t, 3)
            call unpack_full(parvec, k, model%h_t(:, :, t))
         end do
      else if (allocated(model%h)) then
         call unpack_full(parvec, k, model%h)
      end if
      if (allocated(model%l)) call unpack_full(parvec, k, model%l)
      if (allocated(model%q_noise)) model%q = matmul(matmul(model%g, model%q_noise), transpose(model%g))
      if (allocated(model%r_noise)) model%r = matmul(matmul(model%h, model%r_noise), transpose(model%h))
      if (allocated(model%v0_noise)) model%v0 = matmul(matmul(model%l, model%v0_noise), transpose(model%l))
      info = 0
   end subroutine marss_unvectorizeparam

   pure subroutine marss_hatyt(model, kf, result, only_kem)
      type(marss_model), intent(in) :: model !! Model and data used to construct conditional expected observations.
      type(marss_kf_result), intent(in) :: kf !! Filter/smoother moments supplying prior, filtered, and smoothed state &
         !! distributions.
      type(marss_hatyt_result), intent(out) :: result !! Conditional observation means, variances, and state cross moments.
      logical, intent(in), optional :: only_kem !! If true, return only the four moments required by EM; defaults to true.
      real(dp), allocatable :: correction(:, :)
      real(dp), allocatable :: r_inverse(:, :)
      real(dp), allocatable :: r_observed(:, :)
      real(dp), allocatable :: r_columns(:, :)
      real(dp) :: at(size(model%a))
      real(dp) :: covariance(size(model%z, 1), size(model%z, 1))
      real(dp) :: covariance_expected(size(model%z, 1), size(model%z, 1))
      real(dp) :: delta_r(size(model%z, 1), size(model%z, 1))
      real(dp) :: dz(size(model%z, 1), size(model%b, 1))
      real(dp) :: i_missing(size(model%z, 1), size(model%z, 1))
      real(dp) :: mean_y(size(model%z, 1))
      real(dp) :: mean_y_filt(size(model%z, 1))
      real(dp) :: nan_value
      real(dp) :: r_tolerance
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      real(dp) :: v0eff(size(model%v0, 1), size(model%v0, 2))
      real(dp) :: x_previous(size(model%b, 1))
      real(dp) :: y_zero(size(model%z, 1))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      logical :: is_diagonal
      logical :: kem_only
      logical :: missing(size(model%z, 1))
      integer, allocatable :: observed_noisy(:)
      integer :: i
      integer :: info
      integer :: j
      integer :: k
      integer :: m
      integer :: n
      integer :: rank
      integer :: t
      integer :: tt

      result%ok = .false.
      result%info = 0
      kem_only = .true.
      if (present(only_kem)) kem_only = only_kem
      n = size(model%z, 1)
      m = size(model%b, 1)
      tt = size(model%y, 2)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(result%yt(n, tt), result%ot(n, n, tt), result%yxt(n, m, tt))
      allocate(result%yxtp(n, m, tt))
      result%yxtp = nan_value
      if (.not. kem_only) then
         if (.not. allocated(kf%x_pred) .or. .not. allocated(kf%p_pred) .or. &
             .not. allocated(kf%x_filt) .or. .not. allocated(kf%p_filt) .or. &
             .not. allocated(kf%p_lag) .or. .not. allocated(kf%x0_smooth)) then
            result%info = 1
            return
         end if
         allocate(result%var_yt(n, n, tt), result%var_expected_yt(n, n, tt))
         allocate(result%yxt_prev_smooth(n, m, tt))
         allocate(result%ytt1(n, tt), result%ott1(n, n, tt))
         allocate(result%var_ytt1(n, n, tt), result%var_expected_ytt1(n, n, tt))
         allocate(result%yxtt1(n, m, tt), result%ytt(n, tt), result%ott(n, n, tt), result%yxtt(n, m, tt))
         result%var_yt = 0.0_dp
         result%var_expected_yt = 0.0_dp
      end if
      v0eff = marss_v0_effective(model)

      do t = 1, tt
         zt = marss_z_at(model, t)
         at = marss_a_at(model, t)
         rt = marss_r_at(model, t)
         missing = ieee_is_nan(model%y(:, t))
         if (.not. kem_only) then
            result%ytt1(:, t) = matmul(zt, kf%x_pred(:, t)) + at
            result%var_expected_ytt1(:, :, t) = matmul(matmul(zt, kf%p_pred(:, :, t)), transpose(zt))
            result%var_ytt1(:, :, t) = rt + result%var_expected_ytt1(:, :, t)
            result%ott1(:, :, t) = result%var_ytt1(:, :, t) + &
               outer_product(result%ytt1(:, t), result%ytt1(:, t))
            result%yxtt1(:, :, t) = outer_product(result%ytt1(:, t), kf%x_pred(:, t)) + &
               matmul(zt, kf%p_pred(:, :, t))
            if (t == 1) then
               x_previous = kf%x0_smooth
               do i = 1, m
                  if (abs(v0eff(i, i)) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(v0eff(i, i)))) then
                     x_previous(i) = model%x0(i)
                  end if
               end do
            else
               x_previous = kf%x_smooth(:, t - 1)
            end if
         end if

         if (.not. any(missing)) then
            result%yt(:, t) = model%y(:, t)
            result%ot(:, :, t) = outer_product(result%yt(:, t), result%yt(:, t))
            result%yxt(:, :, t) = outer_product(result%yt(:, t), kf%x_smooth(:, t))
            if (t < tt) result%yxtp(:, :, t) = outer_product(result%yt(:, t), kf%x_smooth(:, t + 1))
            if (.not. kem_only) then
               result%ytt(:, t) = result%yt(:, t)
               result%ott(:, :, t) = result%ot(:, :, t)
               result%yxt_prev_smooth(:, :, t) = outer_product(result%yt(:, t), x_previous)
               result%yxtt(:, :, t) = outer_product(result%yt(:, t), kf%x_filt(:, t))
            end if
            cycle
         end if

         y_zero = 0.0_dp
         do i = 1, n
            if (.not. missing(i)) y_zero(i) = model%y(i, t)
         end do
         delta_r = 0.0_dp
         i_missing = 0.0_dp
         do i = 1, n
            delta_r(i, i) = 1.0_dp
            if (missing(i)) i_missing(i, i) = 1.0_dp
         end do

         r_tolerance = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(rt)))
         is_diagonal = .true.
         do j = 1, n
            do i = 1, n
               if (i /= j .and. abs(rt(i, j)) > r_tolerance) is_diagonal = .false.
            end do
         end do
         if (is_diagonal) then
            do i = 1, n
               if (.not. missing(i) .and. abs(rt(i, i)) > r_tolerance) delta_r(i, i) = 0.0_dp
            end do
         else
            k = count((.not. missing) .and. diagonal_nonzero(rt, r_tolerance))
            allocate(observed_noisy(k))
            k = 0
            do i = 1, n
               if (.not. missing(i) .and. abs(rt(i, i)) > r_tolerance) then
                  k = k + 1
                  observed_noisy(k) = i
               end if
            end do
            if (k > 0) then
               allocate(r_observed(k, k), r_columns(n, k))
               do j = 1, k
                  r_columns(:, j) = rt(:, observed_noisy(j))
                  do i = 1, k
                     r_observed(i, j) = rt(observed_noisy(i), observed_noisy(j))
                  end do
               end do
               call psd_inverse(r_observed, r_inverse, rank, info)
               if (info /= 0) then
                  result%info = 100 + t
                  return
               end if
               correction = matmul(r_columns, r_inverse)
               do j = 1, k
                  delta_r(:, observed_noisy(j)) = delta_r(:, observed_noisy(j)) - correction(:, j)
               end do
               deallocate(correction, r_inverse, r_observed, r_columns)
            end if
            deallocate(observed_noisy)
         end if

         mean_y = matmul(zt, kf%x_smooth(:, t)) + at
         result%yt(:, t) = y_zero - matmul(delta_r, y_zero - mean_y)
         dz = matmul(delta_r, zt)
         covariance_expected = matmul(matmul(dz, kf%p_smooth(:, :, t)), transpose(dz))
         covariance = matmul(delta_r, rt) + covariance_expected
         covariance = matmul(matmul(i_missing, covariance), i_missing)
         covariance = 0.5_dp * (covariance + transpose(covariance))
         result%ot(:, :, t) = outer_product(result%yt(:, t), result%yt(:, t)) + covariance
         result%yxt(:, :, t) = outer_product(result%yt(:, t), kf%x_smooth(:, t)) + &
            matmul(dz, kf%p_smooth(:, :, t))
         if (t < tt) then
            result%yxtp(:, :, t) = outer_product(result%yt(:, t), kf%x_smooth(:, t + 1)) + &
               matmul(dz, transpose(kf%p_lag(:, :, t + 1)))
         end if

         if (.not. kem_only) then
            mean_y_filt = matmul(zt, kf%x_filt(:, t)) + at
            result%ytt(:, t) = y_zero - matmul(delta_r, y_zero - mean_y_filt)
            covariance_expected = matmul(matmul(dz, kf%p_filt(:, :, t)), transpose(dz))
            covariance = matmul(delta_r, rt) + covariance_expected
            covariance = matmul(matmul(i_missing, covariance), i_missing)
            covariance = 0.5_dp * (covariance + transpose(covariance))
            result%ott(:, :, t) = outer_product(result%ytt(:, t), result%ytt(:, t)) + covariance
            result%yxt_prev_smooth(:, :, t) = outer_product(result%yt(:, t), x_previous) + &
               matmul(dz, kf%p_lag(:, :, t))
            result%yxtt(:, :, t) = outer_product(result%ytt(:, t), kf%x_filt(:, t)) + &
               matmul(dz, kf%p_filt(:, :, t))
            result%var_expected_yt(:, :, t) = matmul(matmul(i_missing, &
               matmul(matmul(dz, kf%p_smooth(:, :, t)), transpose(dz))), i_missing)
            result%var_yt(:, :, t) = matmul(matmul(i_missing, &
               matmul(delta_r, rt) + matmul(matmul(dz, kf%p_smooth(:, :, t)), transpose(dz))), i_missing)
         end if
      end do
      result%ok = .true.
   end subroutine marss_hatyt

   pure function diagonal_nonzero(a, tolerance) result(nonzero)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose diagonal nonzero pattern is returned.
      real(dp), intent(in) :: tolerance !! Absolute threshold below which diagonal entries are treated as deterministic zeros.
      logical :: nonzero(size(a, 1))
      integer :: i

      do i = 1, size(a, 1)
         nonzero(i) = abs(a(i, i)) > tolerance
      end do
   end function diagonal_nonzero

   pure subroutine marss_residuals(model, kf, result)
      type(marss_model), intent(in) :: model !! Model whose observation residuals are evaluated.
      type(marss_kf_result), intent(in) :: kf !! Filter and smoother output supplying innovations and state estimates.
      type(marss_residual_result), intent(out) :: result !! Innovation, standardized, and smoothed observation residuals.
      real(dp) :: nan_value
      real(dp) :: at(size(model%a))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      integer :: i
      integer :: n
      integer :: t
      integer :: tt

      n = size(model%y, 1)
      tt = size(model%y, 2)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(result%innovations(n, tt), result%standardized_innovations(n, tt), result%smoothed(n, tt))
      result%innovations = kf%innov
      do t = 1, tt
         zt = marss_z_at(model, t)
         at = marss_a_at(model, t)
         do i = 1, n
            if (ieee_is_nan(model%y(i, t))) then
               result%innovations(i, t) = nan_value
               result%standardized_innovations(i, t) = nan_value
               result%smoothed(i, t) = nan_value
            else
               if (kf%sigma(i, i, t) > 0.0_dp) then
                  result%standardized_innovations(i, t) = kf%innov(i, t) / sqrt(kf%sigma(i, i, t))
               else
                  result%standardized_innovations(i, t) = nan_value
               end if
               result%smoothed(i, t) = model%y(i, t) - &
                  dot_product(zt(i, :), kf%x_smooth(:, t)) - at(i)
            end if
         end do
      end do
   end subroutine marss_residuals

   pure subroutine marss_hessian(model, hessian, rel_step, info)
      type(marss_model), intent(in) :: model !! MARSS model at which the numerical log-likelihood Hessian is evaluated.
      real(dp), allocatable, intent(out) :: hessian(:, :) !! Central finite-difference Hessian of the packed parameters.
      real(dp), intent(in), optional :: rel_step !! Relative finite-difference step, defaults to 1e-4.
      integer, intent(out) :: info !! Zero on success, nonzero if a perturbed model cannot be filtered.
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: work(:)
      real(dp) :: f0
      real(dp) :: fmm
      real(dp) :: fmp
      real(dp) :: fpm
      real(dp) :: fpp
      real(dp) :: fm
      real(dp) :: fp
      real(dp) :: hi
      real(dp) :: hj
      real(dp) :: step
      integer :: i
      integer :: j
      integer :: p

      call marss_vectorizeparam(model, theta)
      p = size(theta)
      allocate(hessian(p, p), work(p))
      step = 1.0e-4_dp
      if (present(rel_step)) step = rel_step
      call loglik_at(model, theta, f0, info)
      if (info /= 0) return
      hessian = 0.0_dp
      do i = 1, p
         hi = step * max(1.0_dp, abs(theta(i)))
         work = theta
         work(i) = theta(i) + hi
         call loglik_at(model, work, fp, info)
         if (info /= 0) return
         work(i) = theta(i) - hi
         call loglik_at(model, work, fm, info)
         if (info /= 0) return
         hessian(i, i) = (fp - 2.0_dp*f0 + fm) / (hi*hi)
         do j = i + 1, p
            hj = step * max(1.0_dp, abs(theta(j)))
            work = theta
            work(i) = theta(i) + hi
            work(j) = theta(j) + hj
            call loglik_at(model, work, fpp, info)
            if (info /= 0) return
            work(j) = theta(j) - hj
            call loglik_at(model, work, fpm, info)
            if (info /= 0) return
            work(i) = theta(i) - hi
            work(j) = theta(j) + hj
            call loglik_at(model, work, fmp, info)
            if (info /= 0) return
            work(j) = theta(j) - hj
            call loglik_at(model, work, fmm, info)
            if (info /= 0) return
            hessian(i, j) = (fpp - fpm - fmp + fmm) / (4.0_dp*hi*hj)
            hessian(j, i) = hessian(i, j)
         end do
      end do
      info = 0
   end subroutine marss_hessian

   pure subroutine marss_fisher_i(model, fisher, rel_step, info)
      type(marss_model), intent(in) :: model !! MARSS model at which observed Fisher information is requested.
      real(dp), allocatable, intent(out) :: fisher(:, :) !! Negative numerical Hessian of the log likelihood.
      real(dp), intent(in), optional :: rel_step !! Relative Hessian step size, defaults to 1e-4.
      integer, intent(out) :: info !! Zero on success, nonzero if numerical Hessian evaluation fails.
      real(dp), allocatable :: hessian(:, :)

      if (present(rel_step)) then
         call marss_hessian(model, hessian, rel_step, info)
      else
         call marss_hessian(model, hessian, info=info)
      end if
      if (info == 0) then
         allocate(fisher(size(hessian, 1), size(hessian, 2)))
         fisher = -hessian
      else
         allocate(fisher(0, 0))
      end if
   end subroutine marss_fisher_i

   pure subroutine marss_param_cis(model, alpha, estimate, se, lower, upper, info, rel_step)
      type(marss_model), intent(in) :: model !! Fitted MARSS model whose packed parameter intervals are requested.
      real(dp), intent(in) :: alpha !! Two-sided significance level, for example 0.05 for 95 percent intervals.
      real(dp), allocatable, intent(out) :: estimate(:) !! Packed parameter estimates.
      real(dp), allocatable, intent(out) :: se(:) !! Hessian-based asymptotic standard errors.
      real(dp), allocatable, intent(out) :: lower(:) !! Lower normal-approximation confidence limits.
      real(dp), allocatable, intent(out) :: upper(:) !! Upper normal-approximation confidence limits.
      integer, intent(out) :: info !! Zero on success, nonzero if Fisher information cannot be inverted.
      real(dp), intent(in), optional :: rel_step !! Relative numerical-Hessian step size.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fisher(:, :)
      real(dp) :: zcrit
      integer :: i

      call marss_vectorizeparam(model, estimate)
      if (present(rel_step)) then
         call marss_fisher_i(model, fisher, rel_step, info)
      else
         call marss_fisher_i(model, fisher, info=info)
      end if
      if (info /= 0) return
      call inverse_matrix(fisher, covariance, info)
      if (info /= 0) return
      allocate(se(size(estimate)), lower(size(estimate)), upper(size(estimate)))
      do i = 1, size(estimate)
         se(i) = sqrt(max(covariance(i, i), 0.0_dp))
      end do
      zcrit = normal_quantile(1.0_dp - 0.5_dp*alpha)
      lower = estimate - zcrit*se
      upper = estimate + zcrit*se
   end subroutine marss_param_cis

   pure subroutine loglik_at(template, theta, value, info)
      type(marss_model), intent(in) :: template !! Template supplying data and dimensions for a perturbed model.
      real(dp), intent(in) :: theta(:) !! Packed numerical parameter vector at which to evaluate log likelihood.
      real(dp), intent(out) :: value !! Kalman innovations log likelihood.
      integer, intent(out) :: info !! Zero on success, nonzero if unpacking or filtering fails.
      type(marss_model) :: model
      type(marss_kf_result) :: kf

      call marss_unvectorizeparam(template, theta, model, info)
      if (info /= 0) return
      call marss_kfss(model, kf, smoother=.false.)
      if (.not. kf%ok) then
         info = 100 + kf%info
         return
      end if
      value = kf%loglik
      info = 0
   end subroutine loglik_at

   pure subroutine pack_full(a, vector, k)
      real(dp), intent(in) :: a(:, :) !! Full matrix appended in Fortran column-major order.
      real(dp), intent(inout) :: vector(:) !! Destination parameter vector.
      integer, intent(inout) :: k !! Current last-used position, advanced for each packed value.
      integer :: i
      integer :: j

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            k = k + 1
            vector(k) = a(i, j)
         end do
      end do
   end subroutine pack_full

   pure subroutine unpack_full(vector, k, a)
      real(dp), intent(in) :: vector(:) !! Packed parameter vector containing a full matrix next.
      integer, intent(inout) :: k !! Current last-used position, advanced while unpacking.
      real(dp), intent(out) :: a(:, :) !! Reconstructed full matrix in Fortran column-major order.
      integer :: i
      integer :: j

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            k = k + 1
            a(i, j) = vector(k)
         end do
      end do
   end subroutine unpack_full

   pure subroutine pack_vector(a, vector, k)
      real(dp), intent(in) :: a(:) !! Vector appended to the parameter vector.
      real(dp), intent(inout) :: vector(:) !! Destination parameter vector.
      integer, intent(inout) :: k !! Current last-used position, advanced by the vector length.

      vector(k + 1:k + size(a)) = a
      k = k + size(a)
   end subroutine pack_vector

   pure subroutine unpack_vector(vector, k, a)
      real(dp), intent(in) :: vector(:) !! Packed parameter vector containing a vector next.
      integer, intent(inout) :: k !! Current last-used position, advanced by the vector length.
      real(dp), intent(out) :: a(:) !! Reconstructed vector.

      a = vector(k + 1:k + size(a))
      k = k + size(a)
   end subroutine unpack_vector

   pure subroutine pack_lower(a, vector, k)
      real(dp), intent(in) :: a(:, :) !! Symmetric matrix whose lower triangle is appended column by column.
      real(dp), intent(inout) :: vector(:) !! Destination parameter vector.
      integer, intent(inout) :: k !! Current last-used position, advanced for each packed value.
      integer :: i
      integer :: j

      do j = 1, size(a, 2)
         do i = j, size(a, 1)
            k = k + 1
            vector(k) = a(i, j)
         end do
      end do
   end subroutine pack_lower

   pure subroutine unpack_lower(vector, k, a)
      real(dp), intent(in) :: vector(:) !! Packed parameter vector containing a symmetric lower triangle next.
      integer, intent(inout) :: k !! Current last-used position, advanced while unpacking.
      real(dp), intent(out) :: a(:, :) !! Reconstructed symmetric matrix.
      integer :: i
      integer :: j

      a = 0.0_dp
      do j = 1, size(a, 2)
         do i = j, size(a, 1)
            k = k + 1
            a(i, j) = vector(k)
            a(j, i) = vector(k)
         end do
      end do
   end subroutine unpack_lower

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector in the rank-one product.
      real(dp), intent(in) :: y(:) !! Right vector in the rank-one product.
      real(dp) :: a(size(x), size(y))

      a = spread(x, 2, size(y)) * spread(y, 1, size(x))
   end function outer_product

end module marss_analysis
