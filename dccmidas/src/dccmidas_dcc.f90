! SPDX-License-Identifier: GPL-3.0-only
! DCC, DCC-MIDAS, asymmetric DCC, and DECO kernels translated from dccmidas.
module dccmidas_dcc
   use r_kinds, only : dp
   use rumidas, only : rumidas_lag_weights, RUMIDAS_BETA_LAG, RUMIDAS_ALMON_LAG
   use roll, only : roll_sum
   use dccmidas_types, only : dcc_matrices, DCCMIDAS_SUCCESS, DCCMIDAS_INVALID_INPUT, DCCMIDAS_LINALG_ERROR
   use dccmidas_matrix, only : identity_matrix, outer_product, sample_covariance_series
   use dccmidas_matrix, only : covariance_to_correlation, logdet_quadratic
   implicit none
   private

   public :: dcc_loglik, dcc_mat_est, a_dcc_loglik, a_dcc_mat_est
   public :: dccmidas_loglik, dccmidas_mat_est, a_dccmidas_loglik, a_dccmidas_mat_est
   public :: deco_loglik, deco_mat_est

contains

   pure integer function effective_kc(k_c) result(value)
      integer, intent(in), optional :: k_c !! R-level burn-in count; when supplied, upstream internally adds one.
      value = 2
      if (present(k_c)) value = k_c + 1
   end function effective_kc

   pure subroutine initialize_identity_path(path)
      real(dp), intent(out) :: path(:, :, :) !! Square matrix path to initialize to the identity at every time.
      integer :: t
      path = 0.0_dp
      do t = 1, size(path, 3)
         path(:, :, t) = identity_matrix(size(path, 1))
      end do
   end subroutine initialize_identity_path

   pure subroutine h_from_sd_corr(dt, r, h)
      real(dp), intent(in) :: dt(:, :) !! Diagonal conditional-standard-deviation matrix.
      real(dp), intent(in) :: r(:, :) !! Conditional correlation matrix.
      real(dp), intent(out) :: h(:, :) !! Conditional covariance matrix `Dt R Dt`.
      h = matmul(dt, matmul(r, dt))
   end subroutine h_from_sd_corr

   pure subroutine dcc_loglik(param, res, ll, status, k_c)
      real(dp), intent(in) :: param(:) !! cDCC parameters `(a, b)`.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      real(dp), allocatable, intent(out) :: ll(:) !! Per-time upstream cDCC likelihood contributions.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      integer, intent(in), optional :: k_c !! Initial observations excluded at the R interface; default behavior corresponds to two.
      real(dp), allocatable :: q(:, :, :), r(:, :, :)
      real(dp) :: s(size(res, 1), size(res, 1)), scaled(size(res, 1))
      real(dp) :: logdet, quad, a, b
      integer :: t, start, stat, k, nt, i_diag

      k = size(res, 1)
      nt = size(res, 2)
      allocate(ll(nt), q(k, k, nt), r(k, k, nt))
      ll = 0.0_dp
      if (size(param) /= 2 .or. nt < 2) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      a = param(1)
      b = param(2)
      call sample_covariance_series(res, s, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call initialize_identity_path(q)
      call initialize_identity_path(r)
      start = effective_kc(k_c) + 1
      do t = start, nt
         scaled = sqrt([(q(i_diag, i_diag, t - 1), i_diag = 1, k)]) * res(:, t - 1)
         q(:, :, t) = (1.0_dp - a - b) * s + a * outer_product(scaled) + b * q(:, :, t - 1)
         call covariance_to_correlation(q(:, :, t), r(:, :, t), stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         call logdet_quadratic(r(:, :, t), res(:, t), logdet, quad, stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         ll(t) = -(logdet + quad)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine dcc_loglik

   pure subroutine dcc_mat_est(est_param, res, dt, matrices, status, k_c)
      real(dp), intent(in) :: est_param(:) !! Estimated cDCC parameters `(a, b)`.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      real(dp), intent(in) :: dt(:, :, :) !! Diagonal conditional-standard-deviation matrices.
      type(dcc_matrices), intent(out) :: matrices !! Conditional covariance and correlation paths.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      integer, intent(in), optional :: k_c !! Initial observations excluded at the R interface.
      real(dp), allocatable :: q(:, :, :)
      real(dp) :: s(size(res, 1), size(res, 1)), scaled(size(res, 1))
      real(dp) :: a, b
      integer :: t, start, stat, k, nt, first_count, i_diag

      k = size(res, 1)
      nt = size(res, 2)
      allocate(matrices%h(k, k, nt), matrices%r(k, k, nt), q(k, k, nt))
      if (size(est_param) /= 2 .or. size(dt, 1) /= k .or. size(dt, 2) /= k .or. size(dt, 3) /= nt) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      a = est_param(1)
      b = est_param(2)
      call sample_covariance_series(res, s, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call initialize_identity_path(q)
      call initialize_identity_path(matrices%r)
      matrices%h = 0.0_dp
      first_count = min(effective_kc(k_c), nt)
      do t = 1, first_count
         call h_from_sd_corr(dt(:, :, t), matrices%r(:, :, t), matrices%h(:, :, t))
      end do
      start = effective_kc(k_c) + 1
      do t = start, nt
         scaled = sqrt([(q(i_diag, i_diag, t - 1), i_diag = 1, k)]) * res(:, t - 1)
         q(:, :, t) = (1.0_dp - a - b) * s + a * outer_product(scaled) + b * q(:, :, t - 1)
         call covariance_to_correlation(q(:, :, t), matrices%r(:, :, t), stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         call h_from_sd_corr(dt(:, :, t), matrices%r(:, :, t), matrices%h(:, :, t))
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine dcc_mat_est

   pure subroutine a_dcc_loglik(param, res, ll, status, k_c)
      real(dp), intent(in) :: param(:) !! Asymmetric-DCC parameters `(a, b, g)`.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      real(dp), allocatable, intent(out) :: ll(:) !! Per-time upstream asymmetric-DCC likelihood contributions.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      integer, intent(in), optional :: k_c !! Initial-observation control using the upstream indexing convention.
      real(dp), allocatable :: q(:, :, :), r(:, :, :), eta(:, :)
      real(dp) :: s(size(res, 1), size(res, 1)), nmat(size(res, 1), size(res, 1))
      real(dp) :: logdet, quad, a, b, g
      integer :: t, start, stat, k, nt

      k = size(res, 1)
      nt = size(res, 2)
      allocate(ll(nt), q(k, k, nt), r(k, k, nt), eta(k, nt))
      ll = 0.0_dp
      if (size(param) /= 3 .or. nt < 2) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      a = param(1)
      b = param(2)
      g = param(3)
      eta = min(res, 0.0_dp)
      call sample_covariance_series(res, s, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call sample_covariance_series(eta, nmat, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call initialize_identity_path(q)
      call initialize_identity_path(r)
      start = effective_kc(k_c)
      do t = start, nt
         q(:, :, t) = (1.0_dp - a - b) * s - g * nmat + a * outer_product(res(:, t - 1)) + &
            b * q(:, :, t - 1) + g * outer_product(eta(:, t - 1))
         call covariance_to_correlation(q(:, :, t), r(:, :, t), stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         call logdet_quadratic(r(:, :, t), res(:, t), logdet, quad, stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         ll(t) = -(logdet + quad)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine a_dcc_loglik

   pure subroutine a_dcc_mat_est(est_param, res, dt, matrices, status, k_c)
      real(dp), intent(in) :: est_param(:) !! Estimated asymmetric-DCC parameters `(a, b, g)`.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      real(dp), intent(in) :: dt(:, :, :) !! Diagonal conditional-standard-deviation matrices.
      type(dcc_matrices), intent(out) :: matrices !! Conditional covariance and correlation paths.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      integer, intent(in), optional :: k_c !! Initial-observation control using the upstream convention.
      real(dp), allocatable :: q(:, :, :), eta(:, :)
      real(dp) :: s(size(res, 1), size(res, 1)), nmat(size(res, 1), size(res, 1))
      real(dp) :: a, b, g
      integer :: t, start, stat, k, nt, first_count

      k = size(res, 1)
      nt = size(res, 2)
      allocate(matrices%h(k, k, nt), matrices%r(k, k, nt), q(k, k, nt), eta(k, nt))
      if (size(est_param) /= 3 .or. size(dt, 1) /= k .or. size(dt, 2) /= k .or. size(dt, 3) /= nt) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      a = est_param(1)
      b = est_param(2)
      g = est_param(3)
      eta = min(res, 0.0_dp)
      call sample_covariance_series(res, s, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call sample_covariance_series(eta, nmat, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call initialize_identity_path(q)
      call initialize_identity_path(matrices%r)
      matrices%h = 0.0_dp
      first_count = min(effective_kc(k_c), nt)
      do t = 1, first_count
         call h_from_sd_corr(dt(:, :, t), matrices%r(:, :, t), matrices%h(:, :, t))
      end do
      start = effective_kc(k_c) + 1
      do t = start, nt
         q(:, :, t) = (1.0_dp - a - b) * s - g * nmat + a * outer_product(res(:, t - 1)) + &
            b * q(:, :, t - 1) + g * outer_product(eta(:, t - 1))
         call covariance_to_correlation(q(:, :, t), matrices%r(:, :, t), stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         call h_from_sd_corr(dt(:, :, t), matrices%r(:, :, t), matrices%h(:, :, t))
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine a_dcc_mat_est

   subroutine build_midas_rbar(res, lag_fun, n_c, k_c, w2, r_bar, status)
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      character(len=*), intent(in) :: lag_fun !! Lag weighting name, `Beta` or `Almon`.
      integer, intent(in) :: n_c !! Number of lagged residual realizations used in each local correlation.
      integer, intent(in) :: k_c !! Number of lagged local correlations entering the MIDAS average.
      real(dp), intent(in) :: w2 !! Second Beta/Almon weighting parameter.
      real(dp), allocatable, intent(out) :: r_bar(:, :, :) !! Long-run correlation path with shape `(assets, assets, time)`.
      integer, intent(out) :: status !! Zero on success or a dccmidas/dependency error code.
      real(dp), allocatable :: c_t(:, :, :), weights(:), rolled(:)
      real(dp) :: prod(size(res, 1), size(res, 1)), vdiag(size(res, 1))
      real(dp) :: window(size(res, 1), n_c + 1)
      integer :: i, j, t, k, nt, lag_id, dep_status

      k = size(res, 1)
      nt = size(res, 2)
      allocate(c_t(k, k, nt), r_bar(k, k, nt))
      c_t = 0.0_dp
      r_bar = 1.0_dp
      if (n_c < 0 .or. k_c < 1 .or. nt <= max(n_c, k_c)) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      do t = n_c + 1, nt
         window = res(:, t - n_c:t)
         do i = 1, k
            vdiag(i) = sum(window(i, :) ** 2)
            if (vdiag(i) <= 0.0_dp) then
               status = DCCMIDAS_LINALG_ERROR
               return
            end if
         end do
         prod = matmul(window, transpose(window))
         do j = 1, k
            do i = 1, k
               c_t(i, j, t) = prod(i, j) / sqrt(vdiag(i) * vdiag(j))
            end do
         end do
      end do
      select case (trim(lag_fun))
      case ('Beta', 'beta', 'BETA')
         lag_id = RUMIDAS_BETA_LAG
      case ('Almon', 'almon', 'ALMON')
         lag_id = RUMIDAS_ALMON_LAG
      case default
         status = DCCMIDAS_INVALID_INPUT
         return
      end select
      weights = rumidas_lag_weights(k_c, w2, lag_id, dep_status)
      if (dep_status /= 0 .or. size(weights) /= k_c + 1) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      do j = 1, k
         do i = 1, k
            call roll_sum(c_t(i, j, :), k_c + 1, rolled, weights=weights)
            r_bar(i, j, :) = rolled
         end do
      end do
      do t = 1, min(k_c, nt)
         r_bar(:, :, t) = identity_matrix(k)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine build_midas_rbar

   subroutine dccmidas_loglik(param, res, lag_fun, n_c, k_c, ll, status)
      real(dp), intent(in) :: param(:) !! DCC-MIDAS parameters `(a, b, w2)`.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      character(len=*), intent(in) :: lag_fun !! MIDAS weighting name, `Beta` or `Almon`.
      integer, intent(in) :: n_c !! Residual window length parameter from the R API.
      integer, intent(in) :: k_c !! Long-run-correlation lag count from the R API.
      real(dp), allocatable, intent(out) :: ll(:) !! Per-time upstream DCC-MIDAS likelihood contributions.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp), allocatable :: q(:, :, :), r(:, :, :), r_bar(:, :, :)
      real(dp) :: logdet, quad, a, b
      integer :: t, stat, k, nt

      k = size(res, 1)
      nt = size(res, 2)
      allocate(ll(nt), q(k, k, nt), r(k, k, nt))
      ll = 0.0_dp
      if (size(param) /= 3) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      a = param(1)
      b = param(2)
      call build_midas_rbar(res, lag_fun, n_c, k_c, param(3), r_bar, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call initialize_identity_path(q)
      call initialize_identity_path(r)
      do t = k_c + 1, nt
         q(:, :, t) = (1.0_dp - a - b) * r_bar(:, :, t) + a * outer_product(res(:, t - 1)) + &
            b * q(:, :, t - 1)
         call covariance_to_correlation(q(:, :, t), r(:, :, t), stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         call logdet_quadratic(r(:, :, t), res(:, t), logdet, quad, stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         ll(t) = -(logdet + quad)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine dccmidas_loglik

   subroutine dccmidas_mat_est(est_param, res, dt, lag_fun, n_c, k_c, matrices, status)
      real(dp), intent(in) :: est_param(:) !! Estimated DCC-MIDAS parameters `(a, b, w2)`.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      real(dp), intent(in) :: dt(:, :, :) !! Diagonal conditional-standard-deviation matrices.
      character(len=*), intent(in) :: lag_fun !! MIDAS weighting name, `Beta` or `Almon`.
      integer, intent(in) :: n_c !! Residual window length parameter from the R API.
      integer, intent(in) :: k_c !! Long-run-correlation lag count from the R API.
      type(dcc_matrices), intent(out) :: matrices !! Conditional covariance/correlation and long-run correlation paths.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp), allocatable :: q(:, :, :)
      real(dp) :: a, b
      integer :: t, stat, k, nt

      k = size(res, 1)
      nt = size(res, 2)
      allocate(matrices%h(k, k, nt), matrices%r(k, k, nt), q(k, k, nt))
      if (size(est_param) /= 3 .or. size(dt, 1) /= k .or. size(dt, 2) /= k .or. size(dt, 3) /= nt) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      a = est_param(1)
      b = est_param(2)
      call build_midas_rbar(res, lag_fun, n_c, k_c, est_param(3), matrices%r_bar, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call initialize_identity_path(q)
      call initialize_identity_path(matrices%r)
      matrices%h = 0.0_dp
      do t = 1, min(k_c, nt)
         call h_from_sd_corr(dt(:, :, t), matrices%r(:, :, t), matrices%h(:, :, t))
      end do
      do t = k_c + 1, nt
         q(:, :, t) = (1.0_dp - a - b) * matrices%r_bar(:, :, t) + a * outer_product(res(:, t - 1)) + &
            b * q(:, :, t - 1)
         call covariance_to_correlation(q(:, :, t), matrices%r(:, :, t), stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         call h_from_sd_corr(dt(:, :, t), matrices%r(:, :, t), matrices%h(:, :, t))
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine dccmidas_mat_est

   subroutine a_dccmidas_loglik(param, res, lag_fun, n_c, k_c, ll, status)
      real(dp), intent(in) :: param(:) !! Asymmetric DCC-MIDAS parameters `(a, b, g, w2)`.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      character(len=*), intent(in) :: lag_fun !! MIDAS weighting name, `Beta` or `Almon`.
      integer, intent(in) :: n_c !! Residual window length parameter from the R API.
      integer, intent(in) :: k_c !! Long-run-correlation lag count from the R API.
      real(dp), allocatable, intent(out) :: ll(:) !! Per-time upstream asymmetric DCC-MIDAS likelihood contributions.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp), allocatable :: q(:, :, :), r(:, :, :), r_bar(:, :, :), eta(:, :)
      real(dp) :: nmat(size(res, 1), size(res, 1)), logdet, quad, a, b, g
      integer :: t, stat, k, nt

      k = size(res, 1)
      nt = size(res, 2)
      allocate(ll(nt), q(k, k, nt), r(k, k, nt), eta(k, nt))
      ll = 0.0_dp
      if (size(param) /= 4 .or. nt < 2) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      a = param(1)
      b = param(2)
      g = param(3)
      eta = min(res, 0.0_dp)
      call build_midas_rbar(res, lag_fun, n_c, k_c, param(4), r_bar, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call sample_covariance_series(eta, nmat, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call initialize_identity_path(q)
      call initialize_identity_path(r)
      do t = 2, nt
         q(:, :, t) = (1.0_dp - a - b) * r_bar(:, :, t) - g * nmat + &
            a * outer_product(res(:, t - 1)) + b * q(:, :, t - 1) + g * outer_product(eta(:, t - 1))
         call covariance_to_correlation(q(:, :, t), r(:, :, t), stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         call logdet_quadratic(r(:, :, t), res(:, t), logdet, quad, stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         ll(t) = -(logdet + quad)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine a_dccmidas_loglik

   subroutine a_dccmidas_mat_est(est_param, res, dt, lag_fun, n_c, k_c, matrices, status)
      real(dp), intent(in) :: est_param(:) !! Estimated asymmetric DCC-MIDAS parameters `(a, b, g, w2)`.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      real(dp), intent(in) :: dt(:, :, :) !! Diagonal conditional-standard-deviation matrices.
      character(len=*), intent(in) :: lag_fun !! MIDAS weighting name, `Beta` or `Almon`.
      integer, intent(in) :: n_c !! Residual window length parameter from the R API.
      integer, intent(in) :: k_c !! Long-run-correlation lag count from the R API.
      type(dcc_matrices), intent(out) :: matrices !! Conditional covariance/correlation and long-run correlation paths.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp), allocatable :: q(:, :, :), eta(:, :)
      real(dp) :: nmat(size(res, 1), size(res, 1)), mean_h(size(res, 1), size(res, 1))
      real(dp) :: a, b, g
      integer :: t, stat, k, nt

      k = size(res, 1)
      nt = size(res, 2)
      allocate(matrices%h(k, k, nt), matrices%r(k, k, nt), q(k, k, nt), eta(k, nt))
      if (size(est_param) /= 4 .or. size(dt, 1) /= k .or. size(dt, 2) /= k .or. size(dt, 3) /= nt) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      a = est_param(1)
      b = est_param(2)
      g = est_param(3)
      eta = min(res, 0.0_dp)
      call build_midas_rbar(res, lag_fun, n_c, k_c, est_param(4), matrices%r_bar, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call sample_covariance_series(eta, nmat, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call initialize_identity_path(q)
      call initialize_identity_path(matrices%r)
      matrices%h = 0.0_dp
      do t = k_c + 1, nt
         q(:, :, t) = (1.0_dp - a - b) * matrices%r_bar(:, :, t) - g * nmat + &
            a * outer_product(res(:, t - 1)) + b * q(:, :, t - 1) + g * outer_product(eta(:, t - 1))
         call covariance_to_correlation(q(:, :, t), matrices%r(:, :, t), stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         call h_from_sd_corr(dt(:, :, t), matrices%r(:, :, t), matrices%h(:, :, t))
      end do
      if (k_c > 0 .and. k_c < nt) then
         mean_h = sum(matrices%h(:, :, k_c + 1:nt), dim=3) / real(nt - k_c, dp)
         do t = 1, k_c
            matrices%h(:, :, t) = mean_h
         end do
      end if
      status = DCCMIDAS_SUCCESS
   end subroutine a_dccmidas_mat_est

   pure subroutine deco_loglik(param, res, ll, status, k_c)
      real(dp), intent(in) :: param(:) !! DECO parameters `(a, b)`.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      real(dp), allocatable, intent(out) :: ll(:) !! Per-time upstream DECO likelihood contributions.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      integer, intent(in), optional :: k_c !! Initial observations excluded at the R interface.
      real(dp), allocatable :: q(:, :, :), r(:, :, :)
      real(dp) :: s(size(res, 1), size(res, 1)), scaled(size(res, 1))
      real(dp) :: rdeco(size(res, 1), size(res, 1)), rho, logdet, quad, a, b
      integer :: t, start, stat, i, j, k, nt, i_diag

      k = size(res, 1)
      nt = size(res, 2)
      allocate(ll(nt), q(k, k, nt), r(k, k, nt))
      ll = 0.0_dp
      if (size(param) /= 2 .or. k < 2) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      a = param(1)
      b = param(2)
      call sample_covariance_series(res, s, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call initialize_identity_path(q)
      call initialize_identity_path(r)
      start = effective_kc(k_c) + 1
      do t = start, nt
         scaled = sqrt([(q(i_diag, i_diag, t - 1), i_diag = 1, k)]) * res(:, t - 1)
         q(:, :, t) = (1.0_dp - a - b) * s + a * outer_product(scaled) + b * q(:, :, t - 1)
         call covariance_to_correlation(q(:, :, t), r(:, :, t), stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         rho = (sum(r(:, :, t)) - real(k, dp)) / real(k * (k - 1), dp)
         rdeco = rho
         do j = 1, k
            do i = 1, k
               if (i == j) rdeco(i, j) = 1.0_dp
            end do
         end do
         call logdet_quadratic(rdeco, res(:, t), logdet, quad, stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         ll(t) = -(logdet + quad)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine deco_loglik

   pure subroutine deco_mat_est(est_param, res, dt, matrices, status, k_c)
      real(dp), intent(in) :: est_param(:) !! Estimated DECO parameters `(a, b)`.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      real(dp), intent(in) :: dt(:, :, :) !! Diagonal conditional-standard-deviation matrices.
      type(dcc_matrices), intent(out) :: matrices !! Conditional covariance and equicorrelation paths.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      integer, intent(in), optional :: k_c !! Initial observations excluded at the R interface.
      real(dp), allocatable :: q(:, :, :), rraw(:, :, :)
      real(dp) :: s(size(res, 1), size(res, 1)), scaled(size(res, 1)), rho, a, b
      integer :: t, start, stat, i, j, k, nt, first_count, i_diag

      k = size(res, 1)
      nt = size(res, 2)
      allocate(matrices%h(k, k, nt), matrices%r(k, k, nt), q(k, k, nt), rraw(k, k, nt))
      if (size(est_param) /= 2 .or. k < 2 .or. size(dt, 1) /= k .or. size(dt, 2) /= k .or. size(dt, 3) /= nt) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      a = est_param(1)
      b = est_param(2)
      call sample_covariance_series(res, s, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      call initialize_identity_path(q)
      call initialize_identity_path(rraw)
      call initialize_identity_path(matrices%r)
      matrices%h = 0.0_dp
      first_count = min(effective_kc(k_c), nt)
      do t = 1, first_count
         call h_from_sd_corr(dt(:, :, t), matrices%r(:, :, t), matrices%h(:, :, t))
      end do
      start = effective_kc(k_c) + 1
      do t = start, nt
         scaled = sqrt([(q(i_diag, i_diag, t - 1), i_diag = 1, k)]) * res(:, t - 1)
         q(:, :, t) = (1.0_dp - a - b) * s + a * outer_product(scaled) + b * q(:, :, t - 1)
         call covariance_to_correlation(q(:, :, t), rraw(:, :, t), stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
         rho = (sum(rraw(:, :, t)) - real(k, dp)) / real(k * (k - 1), dp)
         matrices%r(:, :, t) = rho
         do j = 1, k
            do i = 1, k
               if (i == j) matrices%r(i, j, t) = 1.0_dp
            end do
         end do
         call h_from_sd_corr(dt(:, :, t), matrices%r(:, :, t), matrices%h(:, :, t))
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine deco_mat_est

end module dccmidas_dcc
