! SPDX-License-Identifier: GPL-3.0-only
! High-level BEKK and DCC-family fitting workflows translated from dccmidas.
module dccmidas_fit
   use r_kinds, only : dp
   use r_linalg, only : symmetric_eigen, symmetric_eigenvalues
   use maxlik, only : maxlik_problem, maxlik_control, maxlik_result, initialize_problem
   use maxlik, only : set_inequality_constraints, max_lik
   use rugarch, only : garch_fit_result, fit_garch11, fit_gjrgarch11, fit_egarch11, fit_igarch11, fit_csgarch11
   use rugarch, only : dist_norm, dist_std
   use rumidas, only : garch_midas_spec, rumidas_fit_result, fit_garch_midas
   use rumidas, only : RUMIDAS_GM, RUMIDAS_DAGM, RUMIDAS_NORMAL, RUMIDAS_STUDENT_T
   use rumidas, only : RUMIDAS_BETA_LAG, RUMIDAS_ALMON_LAG
   use dccmidas_types, only : dcc_fit_result, bekk_fit_result, dcc_matrices
   use dccmidas_types, only : DCCMIDAS_SUCCESS, DCCMIDAS_INVALID_INPUT, DCCMIDAS_OPTIMIZATION_ERROR
   use dccmidas_bekk, only : sBEKK_loglik, sBEKK_mat_est, dBEKK_loglik, dBEKK_mat_est
   use dccmidas_dcc, only : dcc_loglik, dcc_mat_est, a_dcc_loglik, a_dcc_mat_est
   use dccmidas_dcc, only : dccmidas_loglik, dccmidas_mat_est, a_dccmidas_loglik, a_dccmidas_mat_est
   use dccmidas_dcc, only : deco_loglik, deco_mat_est
   use dccmidas_matrix, only : sample_covariance_series, inv
   use dccmidas_evaluation, only : qmle_sd
   implicit none
   private

   real(dp), allocatable, save :: corr_res(:, :)
   real(dp), allocatable, save :: bekk_ret(:, :)
   character(len=16), save :: corr_name = ''
   character(len=8), save :: corr_lag_fun = 'Beta'
   character(len=8), save :: bekk_name = ''
   integer, save :: corr_n_c = 0
   integer, save :: corr_k_c = 0
   logical, save :: corr_has_k_c = .false.

   public :: dcc_fit, dcc_fit_second_stage, bekk_fit

contains

   subroutine bekk_fit(ret, model, result, status, start)
      real(dp), intent(in) :: ret(:, :) !! Complete returns with shape `(time, assets)`; two to five assets are supported upstream.
      character(len=*), intent(in) :: model !! BEKK model name, `sBEKK` or `dBEKK`.
      type(bekk_fit_result), intent(out) :: result !! Estimated coefficients, QML errors, covariance path, and diagnostics.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp), intent(in), optional :: start(:) !! Optional deterministic optimizer starting vector in upstream parameter order.
      type(maxlik_problem) :: problem
      type(maxlik_control) :: control
      type(maxlik_result) :: fit
      real(dp), allocatable :: initial(:)
      real(dp), allocatable :: a(:, :)
      real(dp), allocatable :: b(:)
      integer :: k
      integer :: m
      integer :: npar
      integer :: stat

      result = bekk_fit_result()
      k = size(ret, 2)
      if (size(ret, 1) < 2 .or. k < 2 .or. k > 5) then
         status = DCCMIDAS_INVALID_INPUT
         result%status = status
         result%message = 'BEKK requires 2 to 5 assets and at least two observations'
         return
      end if
      m = k * (k + 1) / 2
      select case (trim(model))
      case ('sBEKK')
         npar = m + 2
      case ('dBEKK')
         npar = m + 2 * k
      case default
         status = DCCMIDAS_INVALID_INPUT
         result%status = status
         result%message = 'unknown BEKK model'
         return
      end select
      allocate(initial(npar))
      if (present(start)) then
         if (size(start) /= npar) then
            status = DCCMIDAS_INVALID_INPUT
            result%status = status
            result%message = 'BEKK starting vector has the wrong length'
            return
         end if
         initial = start
      else
         initial = merge(0.1_dp, 0.5_dp, trim(model) == 'sBEKK')
      end if
      bekk_ret = ret
      bekk_name = trim(model)
      call initialize_problem(problem, npar, bekk_objective, nobs=size(ret, 1))
      problem%scores => bekk_scores
      if (trim(model) == 'dBEKK') then
         allocate(a(4, npar), b(4))
         a = 0.0_dp
         a(1, m + 1) = 1.0_dp
         a(2, m + k + 1) = 1.0_dp
         a(3, m + 2) = 1.0_dp
         a(4, m + k + 2) = 1.0_dp
         b = [-0.0001_dp, -0.001_dp, -0.0001_dp, -0.001_dp]
         call set_inequality_constraints(problem, a, b, stat)
         if (stat /= 0) then
            status = DCCMIDAS_INVALID_INPUT
            result%status = status
            result%message = 'failed to build dBEKK constraints'
            return
         end if
      end if
      control = maxlik_control()
      control%iterlim = 1000
      control%final_hessian = .true.
      call max_lik(problem, initial, fit, 'bfgs', control)
      if (.not. allocated(fit%estimate)) then
         status = DCCMIDAS_OPTIMIZATION_ERROR
         result%status = status
         result%message = trim(fit%message)
         return
      end if
      result%coefficients = fit%estimate
      result%loglik = fit%maximum
      result%converged = fit%converged
      result%message = trim(fit%message)
      if (allocated(fit%hessian) .and. allocated(fit%gradient_obs)) then
         call qmle_sd(fit%hessian, fit%gradient_obs, result%standard_errors, stat)
      end if
      if (.not. allocated(result%standard_errors) .and. allocated(fit%std_error)) result%standard_errors = fit%std_error
      if (trim(model) == 'sBEKK') then
         call sBEKK_mat_est(result%coefficients, ret, result%h, stat)
      else
         call dBEKK_mat_est(result%coefficients, ret, result%h, stat)
      end if
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         result%status = status
         return
      end if
      status = merge(DCCMIDAS_SUCCESS, DCCMIDAS_OPTIMIZATION_ERROR, fit%converged)
      result%status = status
   end subroutine bekk_fit

   subroutine bekk_objective(x, value, status)
      real(dp), intent(in) :: x(:) !! Candidate BEKK parameter vector.
      real(dp), intent(out) :: value !! Sum of per-observation BEKK log-likelihood contributions.
      integer, intent(out) :: status !! MaxLik-compatible status, zero when the candidate can be evaluated.
      real(dp), allocatable :: ll(:)
      integer :: stat

      if (trim(bekk_name) == 'sBEKK') then
         call sBEKK_loglik(x, bekk_ret, ll, stat)
      else
         call dBEKK_loglik(x, bekk_ret, ll, stat)
      end if
      if (stat == DCCMIDAS_SUCCESS) then
         value = sum(ll)
         status = 0
      else
         value = -huge(1.0_dp)
         status = 1
      end if
   end subroutine bekk_objective

   subroutine bekk_scores(x, scores, status)
      real(dp), intent(in) :: x(:) !! Candidate BEKK parameter vector.
      real(dp), intent(out) :: scores(:, :) !! Observation-wise finite-difference score matrix `(time, parameters)`.
      integer, intent(out) :: status !! MaxLik-compatible status, zero on successful score evaluation.
      real(dp), allocatable :: plus_ll(:)
      real(dp), allocatable :: minus_ll(:)
      real(dp) :: xp(size(x))
      real(dp) :: xm(size(x))
      real(dp) :: step
      integer :: j
      integer :: stat1
      integer :: stat2

      scores = 0.0_dp
      do j = 1, size(x)
         step = epsilon(1.0_dp) ** (1.0_dp / 3.0_dp) * max(1.0_dp, abs(x(j)))
         xp = x
         xm = x
         xp(j) = xp(j) + step
         xm(j) = xm(j) - step
         if (trim(bekk_name) == 'sBEKK') then
            call sBEKK_loglik(xp, bekk_ret, plus_ll, stat1)
            call sBEKK_loglik(xm, bekk_ret, minus_ll, stat2)
         else
            call dBEKK_loglik(xp, bekk_ret, plus_ll, stat1)
            call dBEKK_loglik(xm, bekk_ret, minus_ll, stat2)
         end if
         if (stat1 /= DCCMIDAS_SUCCESS .or. stat2 /= DCCMIDAS_SUCCESS) then
            status = 1
            return
         end if
         scores(:, j) = (plus_ll - minus_ll) / (2.0_dp * step)
      end do
      status = 0
   end subroutine bekk_scores

   subroutine dcc_fit(ret, univ_model, distribution, corr_model, result, status, lag_fun, n_c, k_c, mv, midas_k)
      real(dp), intent(in) :: ret(:, :) !! Complete return matrix with shape `(time, assets)`.
      character(len=*), intent(in) :: univ_model !! First-stage model: rugarch GARCH family or `GM_*`/`DAGM_*`.
      character(len=*), intent(in) :: distribution !! Innovation distribution; `norm` and `std` are supported.
      character(len=*), intent(in) :: corr_model !! Correlation model: `cDCC`, `aDCC`, `DECO`, `DCCMIDAS`, or `ADCCMIDAS`.
      type(dcc_fit_result), intent(out) :: result !! Fitted two-step model and conditional covariance/correlation paths.
      integer, intent(out) :: status !! Zero on success or a dccmidas/dependency error code.
      character(len=*), intent(in), optional :: lag_fun !! MIDAS lag function, `Beta` by default or `Almon`.
      integer, intent(in), optional :: n_c !! Local-correlation residual window required by DCC-MIDAS models.
      integer, intent(in), optional :: k_c !! Long-run correlation lag count; also accepted by non-MIDAS correlation models.
      real(dp), intent(in), optional :: mv(:, :, :) !! Precomputed rumidas lag matrices `(K+1, time, assets)`.
      integer, intent(in), optional :: midas_k !! Number of macro lags represented by the first dimension of `mv` minus one.
      real(dp), allocatable :: sigma(:, :)
      real(dp), allocatable :: standardized(:, :)
      real(dp), allocatable :: dt(:, :, :)
      character(len=8) :: lag_name
      integer :: t
      integer :: i
      integer :: stat
      integer :: nt
      integer :: k

      result = dcc_fit_result()
      nt = size(ret, 1)
      k = size(ret, 2)
      if (nt < 2 .or. k < 2) then
         status = DCCMIDAS_INVALID_INPUT
         result%status = status
         result%message = 'dcc_fit requires at least two assets and two observations'
         return
      end if
      lag_name = 'Beta'
      if (present(lag_fun)) lag_name = trim(lag_fun)
      call fit_univariate_stage(ret, univ_model, distribution, sigma, standardized, stat, lag_name, mv, midas_k)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         result%status = status
         result%message = 'univariate first-stage fit failed'
         return
      end if
      allocate(dt(k, k, nt))
      dt = 0.0_dp
      do t = 1, nt
         do i = 1, k
            dt(i, i, t) = sigma(t, i)
         end do
      end do
      if (present(n_c) .and. present(k_c)) then
         call dcc_fit_second_stage(standardized, dt, corr_model, result, stat, lag_name, n_c, k_c)
      else if (present(k_c)) then
         call dcc_fit_second_stage(standardized, dt, corr_model, result, stat, lag_name, k_c=k_c)
      else
         call dcc_fit_second_stage(standardized, dt, corr_model, result, stat, lag_name)
      end if
      result%conditional_sd = sigma
      result%standardized_residuals = standardized
      status = stat
      result%status = status
   end subroutine dcc_fit

   subroutine fit_univariate_stage(ret, univ_model, distribution, sigma, standardized, status, lag_fun, mv, midas_k)
      real(dp), intent(in) :: ret(:, :) !! Return matrix `(time, assets)` used for independent univariate fits.
      character(len=*), intent(in) :: univ_model !! Upstream first-stage model name.
      character(len=*), intent(in) :: distribution !! Innovation distribution name, `norm` or `std`.
      real(dp), allocatable, intent(out) :: sigma(:, :) !! Conditional standard deviations `(time, assets)`.
      real(dp), allocatable, intent(out) :: standardized(:, :) !! Standardized residuals `(assets, time)`.
      integer, intent(out) :: status !! Zero on success or a dccmidas/dependency error code.
      character(len=*), intent(in) :: lag_fun !! `Beta` or `Almon` for rumidas first-stage fits.
      real(dp), intent(in), optional :: mv(:, :, :) !! Precomputed rumidas lag matrices `(K+1, time, assets)`.
      integer, intent(in), optional :: midas_k !! Number of macro lags used by rumidas GARCH-MIDAS fits.
      type(garch_fit_result) :: gfit
      type(garch_midas_spec) :: spec
      type(rumidas_fit_result) :: mfit
      integer :: distribution_id
      integer :: i
      integer :: stat
      integer :: nt
      integer :: k
      logical :: use_student

      nt = size(ret, 1)
      k = size(ret, 2)
      allocate(sigma(nt, k), standardized(k, nt))
      sigma = 0.0_dp
      standardized = 0.0_dp
      use_student = trim(distribution) == 'std'
      if (trim(distribution) == 'norm') then
         distribution_id = dist_norm
      else if (use_student) then
         distribution_id = dist_std
      else
         status = DCCMIDAS_INVALID_INPUT
         return
      end if

      if (index(trim(univ_model), 'GM_') == 1 .or. index(trim(univ_model), 'DAGM_') == 1) then
         if (.not. present(mv) .or. .not. present(midas_k)) then
            status = DCCMIDAS_INVALID_INPUT
            return
         end if
         if (size(mv, 1) /= midas_k + 1 .or. size(mv, 2) /= nt .or. size(mv, 3) /= k) then
            status = DCCMIDAS_INVALID_INPUT
            return
         end if
         spec = garch_midas_spec()
         spec%k1 = midas_k
         spec%distribution = merge(RUMIDAS_STUDENT_T, RUMIDAS_NORMAL, use_student)
         spec%lag_function = merge(RUMIDAS_ALMON_LAG, RUMIDAS_BETA_LAG, trim(lag_fun) == 'Almon')
         spec%model = merge(RUMIDAS_DAGM, RUMIDAS_GM, index(trim(univ_model), 'DAGM_') == 1)
         spec%skew = index(trim(univ_model), '_skew') > 0 .and. index(trim(univ_model), '_noskew') == 0
         do i = 1, k
            call fit_garch_midas(spec, ret(:, i), mv(:, :, i), mfit, stat)
            if (stat /= 0 .or. .not. allocated(mfit%conditional)) then
               status = DCCMIDAS_OPTIMIZATION_ERROR
               return
            end if
            sigma(:, i) = sqrt(max(mfit%conditional, tiny(1.0_dp)))
            standardized(i, :) = ret(:, i) / sigma(:, i)
         end do
      else
         do i = 1, k
            select case (trim(univ_model))
            case ('sGARCH')
               gfit = fit_garch11(ret(:, i), cond_dist=distribution_id, fit_mean=.false., fit_shape=use_student)
            case ('gjrGARCH')
               gfit = fit_gjrgarch11(ret(:, i), cond_dist=distribution_id, fit_mean=.false., fit_shape=use_student)
            case ('eGARCH')
               gfit = fit_egarch11(ret(:, i), cond_dist=distribution_id, fit_mean=.false., fit_shape=use_student)
            case ('iGARCH')
               gfit = fit_igarch11(ret(:, i), cond_dist=distribution_id, fit_mean=.false., fit_shape=use_student)
            case ('csGARCH')
               gfit = fit_csgarch11(ret(:, i), cond_dist=distribution_id, fit_mean=.false., fit_shape=use_student)
            case default
               status = DCCMIDAS_INVALID_INPUT
               return
            end select
            if (gfit%status /= 0 .or. .not. allocated(gfit%sigma) .or. .not. allocated(gfit%residuals)) then
               status = DCCMIDAS_OPTIMIZATION_ERROR
               return
            end if
            sigma(:, i) = gfit%sigma
            standardized(i, :) = gfit%residuals / max(gfit%sigma, tiny(1.0_dp))
         end do
      end if
      status = DCCMIDAS_SUCCESS
   end subroutine fit_univariate_stage

   subroutine dcc_fit_second_stage(res, dt, corr_model, result, status, lag_fun, n_c, k_c, start)
      real(dp), intent(in) :: res(:, :) !! Standardized residuals with shape `(assets, time)`.
      real(dp), intent(in) :: dt(:, :, :) !! Diagonal conditional-standard-deviation matrices.
      character(len=*), intent(in) :: corr_model !! Correlation model name from the upstream API.
      type(dcc_fit_result), intent(out) :: result !! Fitted correlation-stage parameters and matrix paths.
      integer, intent(out) :: status !! Zero on success or a dccmidas/optimizer error code.
      character(len=*), intent(in), optional :: lag_fun !! MIDAS lag function, `Beta` by default or `Almon`.
      integer, intent(in), optional :: n_c !! DCC-MIDAS local residual window length.
      integer, intent(in), optional :: k_c !! Correlation lag/burn-in control using upstream conventions.
      real(dp), intent(in), optional :: start(:) !! Optional optimizer starting vector.
      type(maxlik_problem) :: problem
      type(maxlik_control) :: control
      type(maxlik_result) :: fit
      real(dp), allocatable :: initial(:)
      real(dp), allocatable :: a(:, :)
      real(dp), allocatable :: b(:)
      real(dp) :: delta
      character(len=8) :: lag_name
      integer :: npar
      integer :: stat

      result = dcc_fit_result()
      lag_name = 'Beta'
      if (present(lag_fun)) lag_name = trim(lag_fun)
      select case (trim(corr_model))
      case ('cDCC', 'DECO')
         npar = 2
      case ('aDCC')
         npar = 3
      case ('DCCMIDAS')
         npar = 3
         if (.not. present(n_c) .or. .not. present(k_c)) then
            status = DCCMIDAS_INVALID_INPUT
            result%status = status
            return
         end if
      case ('ADCCMIDAS')
         npar = 4
         if (.not. present(n_c) .or. .not. present(k_c)) then
            status = DCCMIDAS_INVALID_INPUT
            result%status = status
            return
         end if
      case default
         status = DCCMIDAS_INVALID_INPUT
         result%status = status
         return
      end select
      corr_res = res
      corr_name = trim(corr_model)
      corr_lag_fun = lag_name
      corr_n_c = 0
      if (present(n_c)) corr_n_c = n_c
      corr_has_k_c = present(k_c)
      corr_k_c = 0
      if (present(k_c)) corr_k_c = k_c
      call initialize_problem(problem, npar, corr_objective, nobs=size(res, 2))
      problem%scores => corr_scores
      call build_corr_constraints(corr_model, lag_name, res, a, b, delta, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         result%status = status
         return
      end if
      call set_inequality_constraints(problem, a, b, stat)
      if (stat /= 0) then
         status = DCCMIDAS_INVALID_INPUT
         result%status = status
         return
      end if
      allocate(initial(npar))
      if (present(start)) then
         if (size(start) /= npar) then
            status = DCCMIDAS_INVALID_INPUT
            result%status = status
            return
         end if
         initial = start
      else
         initial = 0.0_dp
         initial(1) = 0.01_dp
         initial(2) = 0.8_dp
         if (trim(corr_model) == 'aDCC' .or. trim(corr_model) == 'ADCCMIDAS') initial(3) = 0.01_dp
         if (trim(corr_model) == 'DCCMIDAS') initial(3) = merge(-0.1_dp, 2.0_dp, trim(lag_name) == 'Almon')
         if (trim(corr_model) == 'ADCCMIDAS') initial(4) = merge(-0.1_dp, 2.0_dp, trim(lag_name) == 'Almon')
      end if
      control = maxlik_control()
      control%iterlim = 1000
      control%final_hessian = .true.
      call max_lik(problem, initial, fit, 'bfgs', control)
      if (.not. allocated(fit%estimate)) then
         status = DCCMIDAS_OPTIMIZATION_ERROR
         result%status = status
         result%message = trim(fit%message)
         return
      end if
      result%coefficients = fit%estimate
      result%loglik = fit%maximum
      result%converged = fit%converged
      result%message = trim(fit%message)
      if (allocated(fit%hessian) .and. allocated(fit%gradient_obs)) then
         call qmle_sd(fit%hessian, fit%gradient_obs, result%standard_errors, stat)
      end if
      if (.not. allocated(result%standard_errors) .and. allocated(fit%std_error)) result%standard_errors = fit%std_error
      call correlation_matrices(fit%estimate, res, dt, result%matrices, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         result%status = status
         return
      end if
      status = merge(DCCMIDAS_SUCCESS, DCCMIDAS_OPTIMIZATION_ERROR, fit%converged)
      result%status = status
   end subroutine dcc_fit_second_stage

   subroutine build_corr_constraints(model, lag_fun, res, a, b, delta, status)
      character(len=*), intent(in) :: model !! Correlation model whose upstream inequality constraints are required.
      character(len=*), intent(in) :: lag_fun !! MIDAS lag function controlling the weight-parameter inequality.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals used to compute the asymmetric stationarity coefficient.
      real(dp), allocatable, intent(out) :: a(:, :) !! MaxLik inequality matrix with constraints `A x + b >= 0`.
      real(dp), allocatable, intent(out) :: b(:) !! MaxLik inequality offsets corresponding to upstream `ineqB` values.
      real(dp), intent(out) :: delta !! Asymmetric stationarity multiplier used by aDCC/ADCCMIDAS.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      integer :: npar
      integer :: nrow
      integer :: w_index
      integer :: stat

      delta = 0.0_dp
      select case (trim(model))
      case ('cDCC', 'DECO')
         npar = 2
         nrow = 3
      case ('DCCMIDAS')
         npar = 3
         nrow = 4
      case ('aDCC')
         npar = 3
         nrow = 5
         call asymmetry_delta(res, delta, stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
      case ('ADCCMIDAS')
         npar = 4
         nrow = 6
         call asymmetry_delta(res, delta, stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            status = stat
            return
         end if
      case default
         status = DCCMIDAS_INVALID_INPUT
         return
      end select
      allocate(a(nrow, npar), b(nrow))
      a = 0.0_dp
      b = 0.0_dp
      a(1, 1) = 1.0_dp
      b(1) = -0.0001_dp
      a(2, 2) = 1.0_dp
      b(2) = -0.001_dp
      if (trim(model) == 'cDCC' .or. trim(model) == 'DECO' .or. trim(model) == 'DCCMIDAS') then
         a(3, 1) = -1.0_dp
         a(3, 2) = -1.0_dp
         b(3) = 0.999_dp
         if (trim(model) == 'DCCMIDAS') then
            w_index = 3
            a(4, w_index) = merge(-1.0_dp, 1.0_dp, trim(lag_fun) == 'Almon')
            b(4) = merge(0.0_dp, -1.001_dp, trim(lag_fun) == 'Almon')
         end if
      else
         a(3, 3) = 1.0_dp
         b(3) = -0.001_dp
         a(4, 1) = -1.0_dp
         a(4, 2) = -1.0_dp
         a(4, 3) = -delta
         b(4) = 0.999_dp
         a(5, 3) = -1.0_dp
         b(5) = merge(0.05_dp, 0.15_dp, trim(model) == 'ADCCMIDAS')
         if (trim(model) == 'ADCCMIDAS') then
            w_index = 4
            a(6, w_index) = merge(-1.0_dp, 1.0_dp, trim(lag_fun) == 'Almon')
            b(6) = merge(0.0_dp, -1.001_dp, trim(lag_fun) == 'Almon')
         end if
      end if
      status = DCCMIDAS_SUCCESS
   end subroutine build_corr_constraints

   subroutine asymmetry_delta(res, delta, status)
      real(dp), intent(in) :: res(:, :) !! Standardized residuals `(assets, time)` used by upstream asymmetric constraints.
      real(dp), intent(out) :: delta !! Largest eigenvalue of the upstream `S^(-1/2) N^(-1) S^(-1/2)` matrix.
      integer, intent(out) :: status !! Zero on success or a dccmidas linear-algebra error code.
      real(dp), allocatable :: eta(:, :)
      real(dp), allocatable :: eval(:)
      real(dp), allocatable :: evec(:, :)
      real(dp), allocatable :: eval_delta(:)
      real(dp), allocatable :: inv_n(:, :)
      real(dp), allocatable :: s_inv_half(:, :)
      real(dp), allocatable :: work(:, :)
      real(dp) :: s(size(res, 1), size(res, 1))
      real(dp) :: nmat(size(res, 1), size(res, 1))
      integer :: i
      integer :: info
      integer :: stat
      integer :: k

      k = size(res, 1)
      allocate(eta(k, size(res, 2)), s_inv_half(k, k), work(k, k))
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
      call symmetric_eigen(s, eval, evec, info)
      if (info /= 0 .or. any(eval <= 0.0_dp)) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      s_inv_half = 0.0_dp
      do i = 1, k
         s_inv_half = s_inv_half + outer_column(evec(:, i), evec(:, i)) / sqrt(eval(i))
      end do
      call inv(nmat, inv_n, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      work = matmul(s_inv_half, matmul(inv_n, s_inv_half))
      call symmetric_eigenvalues(work, eval_delta, info, descending=.true.)
      if (info /= 0 .or. size(eval_delta) < 1) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      delta = eval_delta(1)
      status = DCCMIDAS_SUCCESS
   end subroutine asymmetry_delta

   pure function outer_column(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector in the rank-one matrix product.
      real(dp), intent(in) :: y(:) !! Right vector in the rank-one matrix product.
      real(dp) :: a(size(x), size(y))
      integer :: i
      integer :: j

      do j = 1, size(y)
         do i = 1, size(x)
            a(i, j) = x(i) * y(j)
         end do
      end do
   end function outer_column

   subroutine corr_objective(x, value, status)
      real(dp), intent(in) :: x(:) !! Candidate correlation-stage parameter vector.
      real(dp), intent(out) :: value !! Sum of per-time correlation log-likelihood contributions.
      integer, intent(out) :: status !! MaxLik-compatible status, zero when the candidate can be evaluated.
      real(dp), allocatable :: ll(:)
      integer :: stat

      call corr_ll_vector(x, ll, stat)
      if (stat == DCCMIDAS_SUCCESS) then
         value = sum(ll)
         status = 0
      else
         value = -huge(1.0_dp)
         status = 1
      end if
   end subroutine corr_objective

   subroutine corr_scores(x, scores, status)
      real(dp), intent(in) :: x(:) !! Candidate correlation-stage parameter vector.
      real(dp), intent(out) :: scores(:, :) !! Observation-wise central finite-difference scores `(time, parameters)`.
      integer, intent(out) :: status !! MaxLik-compatible status, zero on successful score evaluation.
      real(dp), allocatable :: plus_ll(:)
      real(dp), allocatable :: minus_ll(:)
      real(dp) :: xp(size(x))
      real(dp) :: xm(size(x))
      real(dp) :: step
      integer :: j
      integer :: stat1
      integer :: stat2

      scores = 0.0_dp
      do j = 1, size(x)
         step = epsilon(1.0_dp) ** (1.0_dp / 3.0_dp) * max(1.0_dp, abs(x(j)))
         xp = x
         xm = x
         xp(j) = xp(j) + step
         xm(j) = xm(j) - step
         call corr_ll_vector(xp, plus_ll, stat1)
         call corr_ll_vector(xm, minus_ll, stat2)
         if (stat1 /= DCCMIDAS_SUCCESS .or. stat2 /= DCCMIDAS_SUCCESS) then
            status = 1
            return
         end if
         scores(:, j) = (plus_ll - minus_ll) / (2.0_dp * step)
      end do
      status = 0
   end subroutine corr_scores

   subroutine corr_ll_vector(x, ll, status)
      real(dp), intent(in) :: x(:) !! Candidate parameters for the active module-level correlation model.
      real(dp), allocatable, intent(out) :: ll(:) !! Per-time log-likelihood contributions for the active model.
      integer, intent(out) :: status !! Zero on success or a dccmidas model-evaluation error code.

      select case (trim(corr_name))
      case ('cDCC')
         if (corr_has_k_c) then
            call dcc_loglik(x, corr_res, ll, status, corr_k_c)
         else
            call dcc_loglik(x, corr_res, ll, status)
         end if
      case ('aDCC')
         if (corr_has_k_c) then
            call a_dcc_loglik(x, corr_res, ll, status, corr_k_c)
         else
            call a_dcc_loglik(x, corr_res, ll, status)
         end if
      case ('DECO')
         if (corr_has_k_c) then
            call deco_loglik(x, corr_res, ll, status, corr_k_c)
         else
            call deco_loglik(x, corr_res, ll, status)
         end if
      case ('DCCMIDAS')
         call dccmidas_loglik(x, corr_res, corr_lag_fun, corr_n_c, corr_k_c, ll, status)
      case ('ADCCMIDAS')
         call a_dccmidas_loglik(x, corr_res, corr_lag_fun, corr_n_c, corr_k_c, ll, status)
      case default
         allocate(ll(size(corr_res, 2)))
         ll = 0.0_dp
         status = DCCMIDAS_INVALID_INPUT
      end select
   end subroutine corr_ll_vector

   subroutine correlation_matrices(x, res, dt, matrices, status)
      real(dp), intent(in) :: x(:) !! Fitted parameter vector for the active correlation model.
      real(dp), intent(in) :: res(:, :) !! Standardized residuals `(assets, time)`.
      real(dp), intent(in) :: dt(:, :, :) !! Diagonal conditional-standard-deviation matrix path.
      type(dcc_matrices), intent(out) :: matrices !! Fitted covariance/correlation and optional long-run-correlation paths.
      integer, intent(out) :: status !! Zero on success or a dccmidas model-evaluation error code.

      select case (trim(corr_name))
      case ('cDCC')
         if (corr_has_k_c) then
            call dcc_mat_est(x, res, dt, matrices, status, corr_k_c)
         else
            call dcc_mat_est(x, res, dt, matrices, status)
         end if
      case ('aDCC')
         if (corr_has_k_c) then
            call a_dcc_mat_est(x, res, dt, matrices, status, corr_k_c)
         else
            call a_dcc_mat_est(x, res, dt, matrices, status)
         end if
      case ('DECO')
         if (corr_has_k_c) then
            call deco_mat_est(x, res, dt, matrices, status, corr_k_c)
         else
            call deco_mat_est(x, res, dt, matrices, status)
         end if
      case ('DCCMIDAS')
         call dccmidas_mat_est(x, res, dt, corr_lag_fun, corr_n_c, corr_k_c, matrices, status)
      case ('ADCCMIDAS')
         call a_dccmidas_mat_est(x, res, dt, corr_lag_fun, corr_n_c, corr_k_c, matrices, status)
      case default
         status = DCCMIDAS_INVALID_INPUT
      end select
   end subroutine correlation_matrices

end module dccmidas_fit
