! SPDX-License-Identifier: GPL-2.0-only
module marss_kfas_bridge
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_kf_result
   use marss_model_ops, only : marss_model_valid
   use marss_parameters, only : marss_b_at, marss_u_at, marss_q_at
   use marss_parameters, only : marss_z_at, marss_a_at, marss_r_at, marss_v0_effective
   use marss_covariance, only : psd_inverse, psd_inverse_logdet
   use kfas, only : kfas_model, kfas_filter_result, kfas_smooth_result
   use kfas, only : kfas_gaussian_filter, kfas_gaussian_loglik, kfas_gaussian_smooth
   implicit none
   private
   public :: marss_kfas
   public :: marss_kfas_loglik
   public :: marss_to_kfas_model

contains

   subroutine marss_kfas(model, result, smoother, return_lag_one)
      type(marss_model), intent(in) :: model !! MARSS model translated to the sibling KFAS exact-diffuse representation.
      type(marss_kf_result), intent(out) :: result !! MARSS-oriented filter, smoother, lag covariance, and likelihood output.
      logical, intent(in), optional :: smoother !! If false, return filtering quantities only; defaults to true.
      logical, intent(in), optional :: return_lag_one !! If false, use the unstacked KFAS state and omit lag-one smoothing.
      type(kfas_model) :: kmodel
      type(kfas_filter_result) :: kfilter
      type(kfas_smooth_result) :: ksmooth
      logical :: do_lag
      logical :: do_smoother
      integer :: info

      result%ok = .false.
      result%info = 0
      result%loglik = -huge(1.0_dp)
      if (.not. marss_model_valid(model)) then
         result%info = 1
         return
      end if
      do_lag = .true.
      if (present(return_lag_one)) do_lag = return_lag_one
      call marss_to_kfas_model(model, kmodel, info, return_lag_one=do_lag)
      if (info /= 0) then
         result%info = 100 + info
         return
      end if
      if (model%diffuse) then
         call kfas_gaussian_filter(kmodel, kfilter, filter_signal=.true., info=info)
      else
         call kfas_gaussian_filter(kmodel, kfilter, tol=0.0_dp, filter_signal=.true., info=info)
      end if
      if (info /= 0) then
         result%info = 200 + info
         return
      end if
      do_smoother = .true.
      if (present(smoother)) do_smoother = smoother
      if (do_smoother) then
         if (model%diffuse) then
            call kfas_gaussian_smooth(kmodel, ksmooth, info=info)
         else
            call kfas_gaussian_smooth(kmodel, ksmooth, tol=0.0_dp, info=info)
         end if
         if (info /= 0) then
            result%info = 300 + info
            return
         end if
         call unpack_kfas_result(model, kfilter, ksmooth, result)
      else
         call unpack_kfas_filter(model, kfilter, result)
      end if
      if (result%info /= 0) return
      result%loglik = kfilter%loglik
      result%ok = .true.
      result%info = 0
   end subroutine marss_kfas

   subroutine marss_kfas_loglik(model, loglik, info, return_lag_one)
      type(marss_model), intent(in) :: model !! MARSS model whose Gaussian likelihood is evaluated by the sibling KFAS backend.
      real(dp), intent(out) :: loglik !! Gaussian log likelihood, or negative huge on model/conversion/KFAS failure.
      integer, intent(out) :: info !! Zero on success; positive values distinguish model, conversion, and KFAS failures.
      logical, intent(in), optional :: return_lag_one !! If true, use the stacked lag-one representation; defaults to false.
      type(kfas_model) :: kmodel
      logical :: do_lag
      integer :: local_info

      loglik = -huge(1.0_dp)
      info = 0
      if (.not. marss_model_valid(model)) then
         info = 1
         return
      end if
      do_lag = .false.
      if (present(return_lag_one)) do_lag = return_lag_one
      call marss_to_kfas_model(model, kmodel, local_info, return_lag_one=do_lag)
      if (local_info /= 0) then
         info = 100 + local_info
         return
      end if
      if (model%diffuse) then
         loglik = kfas_gaussian_loglik(kmodel, info=local_info)
      else
         loglik = kfas_gaussian_loglik(kmodel, tol=0.0_dp, info=local_info)
      end if
      if (local_info /= 0) then
         info = 200 + local_info
         loglik = -huge(1.0_dp)
         return
      end if
   end subroutine marss_kfas_loglik

   subroutine marss_to_kfas_model(model, kmodel, info, return_lag_one)
      type(marss_model), intent(in) :: model !! Valid MARSS model to express as a stacked KFAS Gaussian state-space model.
      type(kfas_model), intent(out) :: kmodel !! KFAS model with intercept and previous-state blocks for lag smoothing.
      integer, intent(out) :: info !! Zero on success; nonzero when model dimensions or diffuse covariance rank are invalid.
      logical, intent(in), optional :: return_lag_one !! If false, build the unstacked KFAS state without lag-one covariance blocks.
      real(dp), allocatable :: inverse(:, :)
      real(dp), allocatable :: gwork(:, :)
      real(dp), allocatable :: qwork(:, :)
      real(dp) :: bt(size(model%b, 1), size(model%b, 2))
      real(dp) :: qt(size(model%q, 1), size(model%q, 2))
      real(dp) :: v0eff(size(model%v0, 1), size(model%v0, 2))
      real(dp) :: x1(size(model%x0))
      real(dp) :: v1(size(model%v0, 1), size(model%v0, 2))
      logical :: do_lag
      logical :: h_varying
      logical :: q_varying
      logical :: r_varying
      logical :: t_varying
      logical :: z_varying
      integer :: disturbance_dim
      integer :: g1
      integer :: i
      integer :: m
      integer :: n
      integer :: nh
      integer :: nq
      integer :: nr
      integer :: nt
      integer :: nz
      integer :: rank
      integer :: s
      integer :: state_dim
      integer :: t
      integer :: tt

      info = 0
      if (.not. marss_model_valid(model)) then
         info = 1
         return
      end if
      n = size(model%y, 1)
      m = size(model%b, 1)
      tt = size(model%y, 2)
      do_lag = .true.
      if (present(return_lag_one)) do_lag = return_lag_one
      state_dim = merge(2 * (m + 1), m + 1, do_lag)
      v0eff = marss_v0_effective(model)
      g1 = process_noise_dimension(model)
      disturbance_dim = merge(2 * g1, g1, do_lag)
      allocate(gwork(m, g1), qwork(g1, g1))

      z_varying = observation_design_is_time_varying(model)
      h_varying = observation_covariance_is_time_varying(model)
      t_varying = transition_is_time_varying(model)
      r_varying = process_loading_is_time_varying(model)
      q_varying = process_covariance_is_time_varying(model)
      kmodel%time_varying = 0
      kmodel%time_varying(1) = merge(1, 0, z_varying)
      kmodel%time_varying(2) = merge(1, 0, h_varying)
      kmodel%time_varying(3) = merge(1, 0, t_varying)
      kmodel%time_varying(4) = merge(1, 0, r_varying)
      kmodel%time_varying(5) = merge(1, 0, q_varying)
      nz = 1 + (tt - 1) * kmodel%time_varying(1)
      nh = 1 + (tt - 1) * kmodel%time_varying(2)
      nt = 1 + (tt - 1) * kmodel%time_varying(3)
      nr = 1 + (tt - 1) * kmodel%time_varying(4)
      nq = 1 + (tt - 1) * kmodel%time_varying(5)

      allocate(kmodel%y(tt, n), kmodel%missing(tt, n))
      allocate(kmodel%z(n, state_dim, nz), kmodel%h(n, n, nh))
      allocate(kmodel%tmat(state_dim, state_dim, nt))
      allocate(kmodel%rmat(state_dim, disturbance_dim, nr), kmodel%q(disturbance_dim, disturbance_dim, nq))
      allocate(kmodel%a1(state_dim))
      allocate(kmodel%p1(state_dim, state_dim))
      allocate(kmodel%p1inf(state_dim, state_dim))
      kmodel%y = 0.0_dp
      kmodel%missing = 0
      kmodel%z = 0.0_dp
      kmodel%h = 0.0_dp
      kmodel%tmat = 0.0_dp
      kmodel%rmat = 0.0_dp
      kmodel%q = 0.0_dp
      kmodel%a1 = 0.0_dp
      kmodel%p1 = 0.0_dp
      kmodel%p1inf = 0.0_dp

      do t = 1, tt
         do i = 1, n
            if (ieee_is_nan(model%y(i, t))) then
               kmodel%missing(t, i) = 1
               kmodel%y(t, i) = 0.0_dp
            else
               kmodel%y(t, i) = model%y(i, t)
            end if
         end do
      end do

      do s = 1, nz
         t = merge(s, 1, z_varying)
         kmodel%z(:, 1:m, s) = marss_z_at(model, t)
         kmodel%z(:, m + 1, s) = marss_a_at(model, t)
      end do
      do s = 1, nh
         t = merge(s, 1, h_varying)
         kmodel%h(:, :, s) = marss_r_at(model, t)
      end do

      if (t_varying) then
         do s = 1, tt - 1
            t = s + 1
            kmodel%tmat(1:m, 1:m, s) = marss_b_at(model, t)
            kmodel%tmat(1:m, m + 1, s) = marss_u_at(model, t)
            kmodel%tmat(m + 1, m + 1, s) = 1.0_dp
            if (do_lag) kmodel%tmat(m + 2:state_dim, 1:m + 1, s) = identity_matrix(m + 1)
         end do
      else
         kmodel%tmat(1:m, 1:m, 1) = marss_b_at(model, 1)
         kmodel%tmat(1:m, m + 1, 1) = marss_u_at(model, 1)
         kmodel%tmat(m + 1, m + 1, 1) = 1.0_dp
         if (do_lag) kmodel%tmat(m + 2:state_dim, 1:m + 1, 1) = identity_matrix(m + 1)
      end if

      if (r_varying) then
         do s = 1, tt - 1
            t = s + 1
            call process_loading_at(model, t, gwork)
            kmodel%rmat(1:m, 1:g1, s) = gwork
         end do
      else
         call process_loading_at(model, 1, gwork)
         kmodel%rmat(1:m, 1:g1, 1) = gwork
      end if

      if (q_varying) then
         do s = 1, tt - 1
            t = s + 1
            call process_covariance_at(model, t, qwork)
            kmodel%q(1:g1, 1:g1, s) = qwork
         end do
      else
         call process_covariance_at(model, 1, qwork)
         kmodel%q(1:g1, 1:g1, 1) = qwork
      end if

      bt = marss_b_at(model, 1)
      qt = marss_q_at(model, 1)
      if (model%tinitx == 0) then
         x1 = matmul(bt, model%x0) + marss_u_at(model, 1)
         v1 = matmul(matmul(bt, v0eff), transpose(bt)) + qt
      else
         x1 = model%x0
         v1 = v0eff
      end if
      kmodel%a1(1:m) = x1
      kmodel%a1(m + 1) = 1.0_dp
      if (do_lag) then
         kmodel%a1(m + 2:2 * m + 1) = model%x0
         kmodel%a1(2 * m + 2) = 1.0_dp
      end if

      if (model%diffuse) then
         kmodel%p1inf(1:m, 1:m) = v1
         if (do_lag .and. model%tinitx == 0) then
            kmodel%p1inf(m + 2:2 * m + 1, m + 2:2 * m + 1) = v0eff
         end if
         call psd_inverse(kmodel%p1inf, inverse, rank, info)
         if (info /= 0) then
            info = 2
            return
         end if
         kmodel%diffuse_rank = rank
      else
         kmodel%p1(1:m, 1:m) = v1
         if (do_lag .and. model%tinitx == 0) then
            kmodel%p1(m + 2:2 * m + 1, m + 2:2 * m + 1) = v0eff
            kmodel%p1(1:m, m + 2:2 * m + 1) = matmul(bt, v0eff)
            kmodel%p1(m + 2:2 * m + 1, 1:m) = transpose(kmodel%p1(1:m, m + 2:2 * m + 1))
         end if
         kmodel%diffuse_rank = 0
      end if
   end subroutine marss_to_kfas_model

   subroutine unpack_kfas_filter(model, kfilter, result)
      type(marss_model), intent(in) :: model !! Original MARSS model used to reconstruct untransformed innovations and gains.
      type(kfas_filter_result), intent(in) :: kfilter !! KFAS exact-diffuse filtering output for the stacked model.
      type(marss_kf_result), intent(inout) :: result !! MARSS filter result populated from KFAS and observation equations.
      integer :: m
      integer :: tt

      m = size(model%b, 1)
      tt = size(model%y, 2)
      allocate(result%x_pred(m, tt), result%p_pred(m, m, tt))
      allocate(result%p_pred_inf(m, m, tt))
      allocate(result%x_filt(m, tt), result%p_filt(m, m, tt))
      allocate(result%x_smooth(m, tt), result%p_smooth(m, m, tt))
      allocate(result%p_lag(m, m, tt), result%x0_smooth(m), result%v0_smooth(m, m))
      result%x_pred = kfilter%a(1:m, 1:tt)
      result%p_pred = kfilter%p(1:m, 1:m, 1:tt)
      result%p_pred_inf = kfilter%pinf(1:m, 1:m, 1:tt)
      result%x_filt = kfilter%att(1:m, :)
      result%p_filt = kfilter%ptt(1:m, 1:m, :)
      result%x_smooth = result%x_filt
      result%p_smooth = result%p_filt
      result%p_lag = 0.0_dp
      result%lag_one_available = .false.
      result%x0_smooth = model%x0
      result%v0_smooth = marss_v0_effective(model)
      result%diffuse_end = kfilter%diffuse_end
      result%diffuse_component = kfilter%diffuse_component
      result%remaining_diffuse_rank = kfilter%remaining_diffuse_rank
      result%observation_covariance_transformed = kfilter%transformed_h
      call reconstruct_observation_filter(model, result)
      if (result%info /= 0) return
      call zero_tiny_covariances_for_deterministic_observations(model, result)
   end subroutine unpack_kfas_filter

   subroutine unpack_kfas_result(model, kfilter, ksmooth, result)
      type(marss_model), intent(in) :: model !! Original MARSS model defining state dimensions and initialization convention.
      type(kfas_filter_result), intent(in) :: kfilter !! KFAS exact-diffuse filtering output for the stacked model.
      type(kfas_smooth_result), intent(in) :: ksmooth !! KFAS exact-diffuse smoothing output including stacked cross covariance.
      type(marss_kf_result), intent(inout) :: result !! MARSS filter/smoother result populated from the KFAS quantities.
      real(dp) :: nan_value
      integer :: m
      integer :: t
      integer :: tt

      call unpack_kfas_filter(model, kfilter, result)
      if (result%info /= 0) return
      m = size(model%b, 1)
      tt = size(model%y, 2)
      result%x_smooth = ksmooth%state(1:m, :)
      result%p_smooth = ksmooth%state_var(1:m, 1:m, :)
      if (size(ksmooth%state, 1) >= 2 * (m + 1)) then
         do t = 1, tt
            result%p_lag(:, :, t) = ksmooth%state_var(1:m, m + 2:2 * m + 1, t)
         end do
         result%lag_one_available = .true.
         if (model%tinitx == 1) then
            nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
            result%p_lag(:, :, 1) = nan_value
         end if
      else
         result%p_lag = 0.0_dp
         result%lag_one_available = .false.
      end if
      result%diffuse_end = ksmooth%diffuse_end
      result%observation_covariance_transformed = result%observation_covariance_transformed .or. &
         ksmooth%transformed_h
      call zero_tiny_covariances_for_deterministic_observations(model, result)
      if (model%tinitx == 0) then
         call reconstruct_initial_smoother(model, result)
      else
         result%x0_smooth = result%x_smooth(:, 1)
         result%v0_smooth = result%p_smooth(:, :, 1)
      end if
   end subroutine unpack_kfas_result

   subroutine reconstruct_initial_smoother(model, result)
      type(marss_model), intent(in) :: model !! MARSS model with tinitx=0 whose x0 smoother moments are reconstructed.
      type(marss_kf_result), intent(inout) :: result !! Smoother output receiving Shumway-Stoffer x0 moments.
      real(dp), allocatable :: vinv(:, :)
      real(dp) :: b1(size(model%b, 1), size(model%b, 2))
      real(dp) :: j0(size(model%b, 1), size(model%b, 2))
      real(dp) :: v0eff(size(model%v0, 1), size(model%v0, 2))
      integer :: info
      integer :: rank

      v0eff = marss_v0_effective(model)
      b1 = marss_b_at(model, 1)
      if (maxval(abs(v0eff)) <= epsilon(1.0_dp)) then
         j0 = 0.0_dp
      else
         call psd_inverse(result%p_pred(:, :, 1), vinv, rank, info)
         if (info /= 0) then
            result%info = 500 + info
            return
         end if
         j0 = matmul(matmul(v0eff, transpose(b1)), vinv)
      end if
      result%x0_smooth = model%x0 + matmul(j0, result%x_smooth(:, 1) - result%x_pred(:, 1))
      result%v0_smooth = v0eff + matmul(matmul(j0, &
         result%p_smooth(:, :, 1) - result%p_pred(:, :, 1)), transpose(j0))
      result%v0_smooth = 0.5_dp * (result%v0_smooth + transpose(result%v0_smooth))
   end subroutine reconstruct_initial_smoother

   subroutine reconstruct_observation_filter(model, result)
      type(marss_model), intent(in) :: model !! MARSS observation equations used to restore quantities before KFAS LDL transforms.
      type(marss_kf_result), intent(inout) :: result !! Result receiving original-coordinate innovations and covariance parts.
      real(dp), allocatable :: finfinv(:, :)
      real(dp), allocatable :: finfobs(:, :)
      real(dp), allocatable :: finv(:, :)
      real(dp), allocatable :: fobs(:, :)
      real(dp), allocatable :: robs(:, :)
      real(dp), allocatable :: zobs(:, :)
      integer, allocatable :: observed(:)
      real(dp) :: at(size(model%a))
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      real(dp) :: logdet
      integer :: i
      integer :: info
      integer :: j
      integer :: k
      integer :: m
      integer :: n
      integer :: rank
      integer :: t
      integer :: tt

      n = size(model%y, 1)
      m = size(model%b, 1)
      tt = size(model%y, 2)
      allocate(result%gain(m, n, tt), result%gain_inf(m, n, tt))
      allocate(result%innov(n, tt), result%sigma(n, n, tt), result%sigma_inf(n, n, tt))
      result%gain = 0.0_dp
      result%gain_inf = 0.0_dp
      result%innov = 0.0_dp
      result%sigma = 0.0_dp
      result%sigma_inf = 0.0_dp
      do t = 1, tt
         zt = marss_z_at(model, t)
         at = marss_a_at(model, t)
         rt = marss_r_at(model, t)
         k = count(.not. ieee_is_nan(model%y(:, t)))
         if (k == 0) cycle
         allocate(observed(k), zobs(k, m), robs(k, k))
         j = 0
         do i = 1, n
            if (.not. ieee_is_nan(model%y(i, t))) then
               j = j + 1
               observed(j) = i
               zobs(j, :) = zt(i, :)
            end if
         end do
         do j = 1, k
            do i = 1, k
               robs(i, j) = rt(observed(i), observed(j))
            end do
         end do
         fobs = matmul(matmul(zobs, result%p_pred(:, :, t)), transpose(zobs)) + robs
         fobs = 0.5_dp * (fobs + transpose(fobs))
         finfobs = matmul(matmul(zobs, result%p_pred_inf(:, :, t)), transpose(zobs))
         finfobs = 0.5_dp * (finfobs + transpose(finfobs))
         call psd_inverse_logdet(fobs, finv, logdet, rank, info)
         if (info /= 0) then
            result%info = 400 + t
            return
         end if
         call psd_inverse(finfobs, finfinv, rank, info)
         if (info /= 0) then
            result%info = 450 + t
            return
         end if
         do i = 1, k
            result%innov(observed(i), t) = model%y(observed(i), t) - &
               dot_product(zt(observed(i), :), result%x_pred(:, t)) - at(observed(i))
            result%gain(:, observed(i), t) = matmul(result%p_pred(:, :, t), &
               matmul(transpose(zobs), finv(:, i)))
            result%gain_inf(:, observed(i), t) = matmul(result%p_pred_inf(:, :, t), &
               matmul(transpose(zobs), finfinv(:, i)))
            do j = 1, k
               result%sigma(observed(i), observed(j), t) = fobs(i, j)
               result%sigma_inf(observed(i), observed(j), t) = finfobs(i, j)
            end do
         end do
         deallocate(observed, zobs, robs, fobs, finfobs, finv, finfinv)
      end do
   end subroutine reconstruct_observation_filter

   subroutine zero_tiny_covariances_for_deterministic_observations(model, result)
      type(marss_model), intent(in) :: model !! Model whose effective observation covariance is checked for deterministic rows.
      type(marss_kf_result), intent(inout) :: result !! Filter/smoother covariances receiving upstream-style epsilon cleanup.
      real(dp), parameter :: tol = epsilon(1.0_dp)

      if (.not. has_deterministic_observation(model)) return
      where (abs(result%p_pred) < tol) result%p_pred = 0.0_dp
      where (abs(result%p_filt) < tol) result%p_filt = 0.0_dp
      where (abs(result%p_smooth) < tol) result%p_smooth = 0.0_dp
      where (abs(result%p_lag) < tol) result%p_lag = 0.0_dp
   end subroutine zero_tiny_covariances_for_deterministic_observations

   pure logical function has_deterministic_observation(model) result(has_zero)
      type(marss_model), intent(in) :: model !! Model tested for a zero diagonal entry in any effective observation covariance.
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      integer :: i
      integer :: t

      has_zero = .false.
      do t = 1, size(model%y, 2)
         rt = marss_r_at(model, t)
         do i = 1, size(rt, 1)
            if (abs(rt(i, i)) <= epsilon(1.0_dp)) then
               has_zero = .true.
               return
            end if
         end do
      end do
   end function has_deterministic_observation

   pure integer function process_noise_dimension(model) result(g1)
      type(marss_model), intent(in) :: model !! Model whose process-disturbance coordinate dimension is requested.

      if (allocated(model%g) .and. allocated(model%q_noise)) then
         g1 = size(model%q_noise, 1)
      else
         g1 = size(model%b, 1)
      end if
   end function process_noise_dimension

   pure subroutine process_loading_at(model, t, value)
      type(marss_model), intent(in) :: model !! Model supplying G(t), or an implicit identity loading for direct Q models.
      integer, intent(in) :: t !! One-based MARSS time index whose process loading is requested.
      real(dp), intent(out) :: value(:, :) !! Process loading with shape state_dimension by disturbance_dimension.

      if (allocated(model%g) .and. allocated(model%q_noise)) then
         if (allocated(model%g_t)) then
            value = model%g_t(:, :, t)
         else
            value = model%g
         end if
      else
         value = identity_matrix(size(model%b, 1))
      end if
   end subroutine process_loading_at

   pure subroutine process_covariance_at(model, t, value)
      type(marss_model), intent(in) :: model !! Model supplying Q(t) in process-disturbance rather than effective state coordinates.
      integer, intent(in) :: t !! One-based MARSS time index whose process-disturbance covariance is requested.
      real(dp), intent(out) :: value(:, :) !! Disturbance covariance with square disturbance-coordinate shape.

      if (allocated(model%g) .and. allocated(model%q_noise)) then
         if (allocated(model%q_noise_t)) then
            value = model%q_noise_t(:, :, t)
         else
            value = model%q_noise
         end if
      else if (allocated(model%q_t)) then
         value = model%q_t(:, :, t)
      else
         value = model%q
      end if
   end subroutine process_covariance_at

   pure logical function observation_design_is_time_varying(model) result(is_varying)
      type(marss_model), intent(in) :: model !! Model tested for a time-dependent combined [Z,A+D*d] observation design.

      is_varying = allocated(model%z_t) .or. allocated(model%a_t) .or. &
         allocated(model%obs_covariates) .or. allocated(model%d_coef_t)
   end function observation_design_is_time_varying

   pure logical function observation_covariance_is_time_varying(model) result(is_varying)
      type(marss_model), intent(in) :: model !! Model tested for a time-dependent effective H(t)R(t)H(t)' covariance.

      is_varying = allocated(model%r_t) .or. allocated(model%h_t) .or. allocated(model%r_noise_t)
   end function observation_covariance_is_time_varying

   pure logical function transition_is_time_varying(model) result(is_varying)
      type(marss_model), intent(in) :: model !! Model tested for a time-dependent combined [B,U+C*c] state transition.

      is_varying = allocated(model%b_t) .or. allocated(model%u_t) .or. &
         allocated(model%state_covariates) .or. allocated(model%c_coef_t)
   end function transition_is_time_varying

   pure logical function process_loading_is_time_varying(model) result(is_varying)
      type(marss_model), intent(in) :: model !! Model tested for a time-dependent process-disturbance loading G(t).

      is_varying = allocated(model%g_t)
   end function process_loading_is_time_varying

   pure logical function process_covariance_is_time_varying(model) result(is_varying)
      type(marss_model), intent(in) :: model !! Model tested for a time-dependent process-disturbance covariance Q(t).

      if (allocated(model%g) .and. allocated(model%q_noise)) then
         is_varying = allocated(model%q_noise_t)
      else
         is_varying = allocated(model%q_t)
      end if
   end function process_covariance_is_time_varying

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n !! Requested order of the identity matrix.
      real(dp) :: a(n, n)
      integer :: i

      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

end module marss_kfas_bridge
