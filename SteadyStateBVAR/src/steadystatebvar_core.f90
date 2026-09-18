! SPDX-License-Identifier: GPL-3.0-or-later
! Derived from SteadyStateBVAR 0.2.0, Copyright (c) 2026 Mark Becker.
module steadystatebvar_core
   use r_kinds, only : dp
   use r_distributions, only : r_qnorm
   use r_quantiles, only : r_quantile_type7
   use r_linalg, only : cholesky_factor, inverse_matrix, solve_spd, solve_system
   use r_linalg, only : spd_inverse_logdet
   use mts, only : var_psi_weights
   use steadystatebvar_types
   implicit none
   private

   public :: bvar_create
   public :: setup_model
   public :: set_priors
   public :: restrict_beta
   public :: fit_model
   public :: forecast_model
   public :: conditional_forecast_model
   public :: irf_model
   public :: ppi
   public :: bvar_summary

contains

   pure subroutine bvar_create(data, model, info)
      real(dp), intent(in) :: data(:, :) !! Endogenous observations, with rows as time and columns as variables.
      type(bvar_model), intent(out) :: model !! Newly initialized BVAR object containing a private copy of `data`.
      integer, intent(out) :: info !! Status code; zero indicates success.

      if (size(data, 1) < 2 .or. size(data, 2) < 1) then
         info = ssbvar_invalid_input
         return
      end if
      allocate(model%data(size(data, 1), size(data, 2)))
      model%data = data
      info = ssbvar_success
   end subroutine bvar_create

   pure subroutine setup_model(model, p, deterministic, dummy, info)
      type(bvar_model), intent(inout) :: model !! BVAR object whose lagged design matrices and OLS estimates are prepared.
      integer, intent(in) :: p !! Positive VAR lag order, strictly smaller than the number of observations.
      character(len=*), intent(in), optional :: deterministic !! One of constant, constant_and_dummy, or constant_and_trend.
      real(dp), intent(in), optional :: dummy(:) !! Dummy series with one value per observation for constant_and_dummy models.
      integer, intent(out) :: info !! Status code; zero indicates success.
      character(len=24) :: det
      integer :: n_obs, k, n, q, kp, t, lag, i, df
      real(dp), allocatable :: dt(:, :), z(:, :), xtx(:, :), xty(:, :), beta_hat(:, :)
      real(dp), allocatable :: residuals(:, :), a_l(:, :, :), a_long(:, :), c_hat(:, :)
      real(dp), allocatable :: z_ar(:, :), xtx_ar(:, :), xty_ar(:), beta_ar(:), u_ar(:)

      info = ssbvar_success
      if (.not. allocated(model%data)) then
         info = ssbvar_not_ready
         return
      end if
      n_obs = size(model%data, 1)
      k = size(model%data, 2)
      if (p <= 0 .or. p >= n_obs) then
         info = ssbvar_invalid_input
         return
      end if
      if (present(dummy)) then
         if (size(dummy) /= n_obs) then
            info = ssbvar_invalid_input
            return
         end if
      end if

      det = "constant"
      if (present(deterministic)) det = trim(deterministic)
      select case (trim(det))
      case ("constant")
         q = 1
      case ("constant_and_dummy")
         q = 2
         if (.not. present(dummy)) then
            info = ssbvar_invalid_input
            return
         end if
      case ("constant_and_trend")
         q = 2
      case default
         info = ssbvar_invalid_input
         return
      end select

      n = n_obs - p
      kp = k * p
      df = n - kp - q
      if (df <= 0) then
         info = ssbvar_invalid_input
         return
      end if

      model%has_setup = .false.
      model%has_priors = .false.
      model%has_fit = .false.
      model%setup = bvar_setup()
      model%priors = bvar_priors()
      model%fit = bvar_fit()
      allocate(dt(n_obs, q))
      dt(:, 1) = 1.0_dp
      if (q == 2) then
         if (trim(det) == "constant_and_dummy") then
            dt(:, 2) = dummy
         else
            do t = 1, n_obs
               dt(t, 2) = real(t, dp)
            end do
         end if
      end if

      model%setup%n = n
      model%setup%k = k
      model%setup%p = p
      model%setup%q = q
      model%setup%n_free_params_a = k * (k - 1) / 2
      model%setup%deterministic = det
      allocate(model%setup%y(n, k), model%setup%x(n, q), model%setup%w(n, kp))
      allocate(model%setup%q_lag(n, q * p), model%setup%d(n, q), model%setup%dt(n_obs, q))
      model%setup%y = model%data(p + 1:n_obs, :)
      model%setup%x = dt(p + 1:n_obs, :)
      model%setup%d = model%setup%x
      model%setup%dt = dt
      do t = 1, n
         do lag = 1, p
            model%setup%w(t, (lag - 1) * k + 1:lag * k) = model%data(p + t - lag, :)
            model%setup%q_lag(t, (lag - 1) * q + 1:lag * q) = dt(p + t - lag, :)
         end do
      end do
      if (present(dummy)) then
         allocate(model%setup%dummy(n_obs))
         model%setup%dummy = dummy
      end if

      allocate(z(n, kp + q))
      z(:, 1:kp) = model%setup%w
      z(:, kp + 1:kp + q) = model%setup%x
      allocate(xtx(kp + q, kp + q), xty(kp + q, k), beta_hat(kp + q, k))
      xtx = matmul(transpose(z), z)
      xty = matmul(transpose(z), model%setup%y)
      call solve_system(xtx, xty, beta_hat, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if

      allocate(residuals(n, k))
      residuals = model%setup%y - matmul(z, beta_hat)
      allocate(model%setup%sigma_u_ols(k, k))
      model%setup%sigma_u_ols = matmul(transpose(residuals), residuals) / real(df, dp)
      allocate(model%setup%beta_ols(kp, k))
      model%setup%beta_ols = beta_hat(1:kp, :)

      allocate(a_l(k, k, p), a_long(k, k), c_hat(k, q))
      do lag = 1, p
         a_l(:, :, lag) = transpose(beta_hat((lag - 1) * k + 1:lag * k, :))
      end do
      a_long = identity_matrix(k)
      do lag = 1, p
         a_long = a_long - a_l(:, :, lag)
      end do
      if (q == 1) then
         c_hat(:, 1) = beta_hat(kp + 1, :)
      else
         c_hat = transpose(beta_hat(kp + 1:kp + q, :))
      end if
      allocate(model%setup%psi_ols(k, q))
      call solve_system(a_long, c_hat, model%setup%psi_ols, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if

      allocate(model%setup%sigma_ar(k, k))
      model%setup%sigma_ar = 0.0_dp
      do i = 1, k
         allocate(z_ar(n, p + q), xtx_ar(p + q, p + q), xty_ar(p + q))
         allocate(beta_ar(p + q), u_ar(n))
         do t = 1, n
            do lag = 1, p
               z_ar(t, lag) = model%data(p + t - lag, i)
            end do
         end do
         z_ar(:, p + 1:p + q) = model%setup%x
         xtx_ar = matmul(transpose(z_ar), z_ar)
         xty_ar = matmul(transpose(z_ar), model%setup%y(:, i))
         call solve_system(xtx_ar, xty_ar, beta_ar, info)
         if (info /= 0) then
            info = ssbvar_linalg_failure
            return
         end if
         u_ar = model%setup%y(:, i) - matmul(z_ar, beta_ar)
         model%setup%sigma_ar(i, i) = sum(u_ar * u_ar) / real(n - p - q, dp)
         deallocate(z_ar, xtx_ar, xty_ar, beta_ar, u_ar)
      end do

      model%has_setup = .true.
      model%has_priors = .false.
      model%has_fit = .false.
      info = ssbvar_success
   end subroutine setup_model

   pure subroutine set_priors(model, lambda_1, lambda_2, lambda_3, first_own_lag_prior_mean, &
      theta_psi, omega_psi, jeffreys, sv, sv_type, theta_a, omega_a, theta_log_lambda_1, &
      omega_log_lambda_1, alpha_phi, beta_phi, theta_gamma_0, omega_gamma_0, theta_gamma_1, &
      omega_gamma_1, m_phi, v_phi, info)
      type(bvar_model), intent(inout) :: model !! Set-up BVAR object receiving Minnesota, steady-state, and optional SV priors.
      real(dp), intent(in), optional :: lambda_1 !! Minnesota overall tightness; default 0.2 and required positive.
      real(dp), intent(in), optional :: lambda_2 !! Minnesota cross-equation tightness; default 0.5 and required positive.
      real(dp), intent(in), optional :: lambda_3 !! Minnesota lag-decay exponent; default 1.0 and required positive.
      real(dp), intent(in), optional :: first_own_lag_prior_mean(:) !! Length-k means for first own-lag coefficients; default zero.
      real(dp), intent(in), optional :: theta_psi(:) !! Prior mean of vec(Psi), length k*q; defaults to vec(Psi_OLS).
      real(dp), intent(in), optional :: omega_psi(:, :) !! Prior covariance of vec(Psi), shape (k*q,k*q); default identity.
      logical, intent(in), optional :: jeffreys !! Use Jeffreys covariance prior when true; ignored and forced false for SV models.
      logical, intent(in), optional :: sv !! Enable stochastic volatility when true; default false.
      character(len=*), intent(in), optional :: sv_type !! Stochastic-volatility type, exactly RW or AR1 when `sv` is true.
      real(dp), intent(in), optional :: theta_a(:) !! Prior means of lower-triangular free A parameters in row-major order.
      real(dp), intent(in), optional :: omega_a(:, :) !! Prior covariance of free A parameters; upstream Stan uses its diagonal.
      real(dp), intent(in), optional :: theta_log_lambda_1(:) !! Length-k prior mean for the initial log-volatility state.
      real(dp), intent(in), optional :: omega_log_lambda_1(:, :) !! Initial log-volatility covariance; Stan uses its diagonal.
      real(dp), intent(in), optional :: alpha_phi(:) !! Positive length-k RW inverse-gamma shapes for volatility innovations.
      real(dp), intent(in), optional :: beta_phi(:) !! Positive length-k RW inverse-gamma scales for volatility innovations.
      real(dp), intent(in), optional :: theta_gamma_0(:) !! AR1 length-k prior means for log-volatility intercepts.
      real(dp), intent(in), optional :: omega_gamma_0(:, :) !! AR1 intercept prior covariance; upstream Stan uses its diagonal.
      real(dp), intent(in), optional :: theta_gamma_1(:) !! AR1 length-k prior means for log-volatility slopes.
      real(dp), intent(in), optional :: omega_gamma_1(:, :) !! AR1 slope prior covariance; upstream Stan uses its diagonal.
      integer, intent(in), optional :: m_phi !! AR1 inverse-Wishart degrees of freedom for Phi; required to be at least k.
      real(dp), intent(in), optional :: v_phi(:, :) !! AR1 inverse-Wishart scale matrix Phi, with shape (k,k).
      integer, intent(out) :: info !! Status code; zero indicates successful prior construction and validation.
      integer :: k, p, q, kp, d_beta, d_psi, n_free, i, j, lag, idx
      real(dp) :: l1, l2, l3, variance
      real(dp), allocatable :: sigma(:), own_mean(:)
      logical :: use_jeffreys, use_sv

      if (.not. model%has_setup) then
         info = ssbvar_not_ready
         return
      end if
      l1 = 0.2_dp
      l2 = 0.5_dp
      l3 = 1.0_dp
      if (present(lambda_1)) l1 = lambda_1
      if (present(lambda_2)) l2 = lambda_2
      if (present(lambda_3)) l3 = lambda_3
      if (l1 <= 0.0_dp .or. l2 <= 0.0_dp .or. l3 <= 0.0_dp) then
         info = ssbvar_invalid_input
         return
      end if

      k = model%setup%k
      p = model%setup%p
      q = model%setup%q
      kp = k * p
      d_beta = kp * k
      d_psi = k * q
      n_free = model%setup%n_free_params_a
      model%has_priors = .false.
      model%has_fit = .false.
      model%priors = bvar_priors()
      model%fit = bvar_fit()
      allocate(sigma(k), own_mean(k))
      do i = 1, k
         if (model%setup%sigma_ar(i, i) <= 0.0_dp) then
            info = ssbvar_invalid_input
            return
         end if
         sigma(i) = sqrt(model%setup%sigma_ar(i, i))
      end do
      own_mean = 0.0_dp
      if (present(first_own_lag_prior_mean)) then
         if (size(first_own_lag_prior_mean) /= k) then
            info = ssbvar_invalid_input
            return
         end if
         own_mean = first_own_lag_prior_mean
      end if

      allocate(model%priors%theta_beta(d_beta), model%priors%omega_beta(d_beta, d_beta))
      model%priors%theta_beta = 0.0_dp
      model%priors%omega_beta = 0.0_dp
      do i = 1, k
         idx = (i - 1) * kp + i
         model%priors%theta_beta(idx) = own_mean(i)
         do lag = 1, p
            do j = 1, k
               idx = (i - 1) * kp + (lag - 1) * k + j
               if (i == j) then
                  variance = (l1 / real(lag, dp) ** l3) ** 2
               else
                  variance = (l1 * l2 * sigma(i) / &
                     (real(lag, dp) ** l3 * sigma(j))) ** 2
               end if
               model%priors%omega_beta(idx, idx) = variance
            end do
         end do
      end do

      allocate(model%priors%theta_psi(d_psi), model%priors%omega_psi(d_psi, d_psi))
      if (present(theta_psi)) then
         if (size(theta_psi) /= d_psi) then
            info = ssbvar_invalid_input
            return
         end if
         model%priors%theta_psi = theta_psi
      else
         model%priors%theta_psi = reshape(model%setup%psi_ols, [d_psi])
      end if
      model%priors%omega_psi = 0.0_dp
      if (present(omega_psi)) then
         if (size(omega_psi, 1) /= d_psi .or. size(omega_psi, 2) /= d_psi) then
            info = ssbvar_invalid_input
            return
         end if
         model%priors%omega_psi = omega_psi
      else
         do i = 1, d_psi
            model%priors%omega_psi(i, i) = 1.0_dp
         end do
      end if
      do i = 1, d_beta
         if (model%priors%omega_beta(i, i) <= 0.0_dp) then
            info = ssbvar_invalid_input
            return
         end if
      end do
      do i = 1, d_psi
         if (model%priors%omega_psi(i, i) <= 0.0_dp) then
            info = ssbvar_invalid_input
            return
         end if
      end do

      use_jeffreys = .true.
      if (present(jeffreys)) use_jeffreys = jeffreys
      use_sv = .false.
      if (present(sv)) use_sv = sv
      if (use_sv) then
         if (k < 2) then
            info = ssbvar_invalid_input
            return
         end if
         if (.not. present(sv_type)) then
            info = ssbvar_invalid_input
            return
         end if
         if (trim(sv_type) /= "RW" .and. trim(sv_type) /= "AR1") then
            info = ssbvar_invalid_input
            return
         end if
         if (.not. present(theta_a) .or. .not. present(omega_a) .or. &
             .not. present(theta_log_lambda_1) .or. .not. present(omega_log_lambda_1)) then
            info = ssbvar_invalid_input
            return
         end if
         if (size(theta_a) /= n_free) then
            info = ssbvar_invalid_input
            return
         end if
         if (size(omega_a, 1) /= n_free .or. size(omega_a, 2) /= n_free) then
            info = ssbvar_invalid_input
            return
         end if
         if (size(theta_log_lambda_1) /= k) then
            info = ssbvar_invalid_input
            return
         end if
         if (size(omega_log_lambda_1, 1) /= k .or. size(omega_log_lambda_1, 2) /= k) then
            info = ssbvar_invalid_input
            return
         end if
         do i = 1, n_free
            if (omega_a(i, i) <= 0.0_dp) then
               info = ssbvar_invalid_input
               return
            end if
         end do
         do i = 1, k
            if (omega_log_lambda_1(i, i) <= 0.0_dp) then
               info = ssbvar_invalid_input
               return
            end if
         end do

         allocate(model%priors%sv_priors%theta_a(n_free))
         allocate(model%priors%sv_priors%omega_a(n_free, n_free))
         allocate(model%priors%sv_priors%theta_log_lambda_1(k))
         allocate(model%priors%sv_priors%omega_log_lambda_1(k, k))
         model%priors%sv_priors%theta_a = theta_a
         model%priors%sv_priors%omega_a = omega_a
         model%priors%sv_priors%theta_log_lambda_1 = theta_log_lambda_1
         model%priors%sv_priors%omega_log_lambda_1 = omega_log_lambda_1

         if (trim(sv_type) == "RW") then
            if (.not. present(alpha_phi) .or. .not. present(beta_phi)) then
               info = ssbvar_invalid_input
               return
            end if
            if (size(alpha_phi) /= k .or. size(beta_phi) /= k) then
               info = ssbvar_invalid_input
               return
            end if
            if (any(alpha_phi <= 0.0_dp) .or. any(beta_phi <= 0.0_dp)) then
               info = ssbvar_invalid_input
               return
            end if
            allocate(model%priors%sv_priors%alpha_phi(k), model%priors%sv_priors%beta_phi(k))
            model%priors%sv_priors%alpha_phi = alpha_phi
            model%priors%sv_priors%beta_phi = beta_phi
         else
            if (.not. present(theta_gamma_0) .or. .not. present(omega_gamma_0) .or. &
                .not. present(theta_gamma_1) .or. .not. present(omega_gamma_1) .or. &
                .not. present(m_phi) .or. .not. present(v_phi)) then
               info = ssbvar_invalid_input
               return
            end if
            if (size(theta_gamma_0) /= k .or. size(theta_gamma_1) /= k) then
               info = ssbvar_invalid_input
               return
            end if
            if (size(omega_gamma_0, 1) /= k .or. size(omega_gamma_0, 2) /= k) then
               info = ssbvar_invalid_input
               return
            end if
            if (size(omega_gamma_1, 1) /= k .or. size(omega_gamma_1, 2) /= k) then
               info = ssbvar_invalid_input
               return
            end if
            if (m_phi < k) then
               info = ssbvar_invalid_input
               return
            end if
            if (size(v_phi, 1) /= k .or. size(v_phi, 2) /= k) then
               info = ssbvar_invalid_input
               return
            end if
            do i = 1, k
               if (omega_gamma_0(i, i) <= 0.0_dp .or. omega_gamma_1(i, i) <= 0.0_dp) then
                  info = ssbvar_invalid_input
                  return
               end if
            end do
            allocate(model%priors%sv_priors%theta_gamma_0(k))
            allocate(model%priors%sv_priors%omega_gamma_0(k, k))
            allocate(model%priors%sv_priors%theta_gamma_1(k))
            allocate(model%priors%sv_priors%omega_gamma_1(k, k))
            allocate(model%priors%sv_priors%v_phi(k, k))
            model%priors%sv_priors%theta_gamma_0 = theta_gamma_0
            model%priors%sv_priors%omega_gamma_0 = omega_gamma_0
            model%priors%sv_priors%theta_gamma_1 = theta_gamma_1
            model%priors%sv_priors%omega_gamma_1 = omega_gamma_1
            model%priors%sv_priors%m_phi = m_phi
            model%priors%sv_priors%v_phi = v_phi
         end if
      end if

      model%priors%jeffreys = use_jeffreys
      model%priors%sv = use_sv
      model%priors%sv_type = ""
      if (present(sv_type)) model%priors%sv_type = trim(sv_type)
      allocate(model%priors%sigma_ar(k, k))
      model%priors%sigma_ar = model%setup%sigma_ar
      if (.not. use_jeffreys .and. .not. use_sv) then
         model%priors%m = k + 2
         allocate(model%priors%v(k, k))
         model%priors%v = real(model%priors%m - k - 1, dp) * model%setup%sigma_u_ols
      end if
      if (use_sv) model%priors%jeffreys = .false.

      model%has_priors = .true.
      model%has_fit = .false.
      info = ssbvar_success
   end subroutine set_priors

   pure subroutine restrict_beta(model, restriction_matrix, info)
      type(bvar_model), intent(inout) :: model !! BVAR object whose beta prior variance is tightened at zero-restriction locations.
      real(dp), intent(in) :: restriction_matrix(:, :) !! Shape (k*p,k); zeros mark coefficients restricted near zero.
      integer, intent(out) :: info !! Status code; zero indicates successful restriction of the prior covariance.
      integer :: kp, k, i, j, idx

      if (.not. model%has_priors) then
         info = ssbvar_not_ready
         return
      end if
      kp = model%setup%k * model%setup%p
      k = model%setup%k
      if (size(restriction_matrix, 1) /= kp .or. size(restriction_matrix, 2) /= k) then
         info = ssbvar_invalid_input
         return
      end if
      if (allocated(model%priors%restriction)) deallocate(model%priors%restriction)
      allocate(model%priors%restriction(kp, k))
      model%priors%restriction = restriction_matrix
      do j = 1, k
         do i = 1, kp
            if (abs(restriction_matrix(i, j)) <= 0.0_dp) then
               idx = (j - 1) * kp + i
               model%priors%omega_beta(idx, idx) = 1.0e-7_dp
            end if
         end do
      end do
      model%has_fit = .false.
      info = ssbvar_success
   end subroutine restrict_beta

   pure subroutine ppi(lower_bound, upper_bound, interval, annualized_growthrate, freq, mean_value, variance, info)
      real(dp), intent(in) :: lower_bound !! Lower bound of the symmetric normal prior probability interval.
      real(dp), intent(in) :: upper_bound !! Upper endpoint of the normal interval; must exceed the lower endpoint.
      real(dp), intent(in), optional :: interval !! Interval probability mass in (0,1); default 0.95.
      logical, intent(in), optional :: annualized_growthrate !! Divide elicited mean and standard deviation by `freq` when true.
      integer, intent(in), optional :: freq !! Positive observations-per-year frequency used for annualized growth-rate priors.
      real(dp), intent(out) :: mean_value !! Implied normal prior mean on the model scale.
      real(dp), intent(out) :: variance !! Implied normal prior variance on the model scale.
      integer, intent(out) :: info !! Status code; zero indicates valid interval conversion.
      real(dp) :: mass, alpha, z, sigma
      logical :: annualized

      mean_value = 0.0_dp
      variance = 0.0_dp
      mass = 0.95_dp
      if (present(interval)) mass = interval
      annualized = .false.
      if (present(annualized_growthrate)) annualized = annualized_growthrate
      if (lower_bound >= upper_bound .or. mass <= 0.0_dp .or. mass >= 1.0_dp) then
         info = ssbvar_invalid_input
         return
      end if
      if (annualized) then
         if (.not. present(freq)) then
            info = ssbvar_invalid_input
            return
         end if
         if (freq <= 0) then
            info = ssbvar_invalid_input
            return
         end if
      end if
      alpha = 1.0_dp - mass
      z = r_qnorm(1.0_dp - alpha / 2.0_dp)
      mean_value = (upper_bound + lower_bound) / 2.0_dp
      sigma = (upper_bound - lower_bound) / (2.0_dp * z)
      if (annualized) then
         mean_value = mean_value / real(freq, dp)
         sigma = sigma / real(freq, dp)
      end if
      variance = sigma * sigma
      info = ssbvar_success
   end subroutine ppi

   subroutine fit_model(model, h, n_draws, burnin, seed, d_pred, sv_state_scale, sv_parameter_scale, info)
      type(bvar_model), intent(inout) :: model !! BVAR object with setup and priors; receives posterior and predictive draws.
      integer, intent(in), optional :: h !! Positive forecast horizon; default one.
      integer, intent(in), optional :: n_draws !! Number of retained posterior draws; default 1000.
      integer, intent(in), optional :: burnin !! Number of discarded MCMC iterations; default 500.
      integer, intent(in), optional :: seed !! Deterministic seed used to initialize the Fortran intrinsic RNG; default 12345.
      real(dp), intent(in), optional :: d_pred(:, :) !! Future deterministic matrix, shape (H,q); default extends observed terms.
      real(dp), intent(in), optional :: sv_state_scale !! Positive latent log-volatility proposal SD; default 0.35.
      real(dp), intent(in), optional :: sv_parameter_scale !! Positive AR1 gamma proposal SD; default 0.05.
      integer, intent(out) :: info !! Status code; zero indicates posterior sampling and predictive simulation succeeded.
      integer :: horizon, draws, warmup, rng_seed, total_iter, iter, kept
      integer :: k, p, q, kp
      real(dp) :: state_scale, parameter_scale
      real(dp), allocatable :: beta(:, :), psi(:, :), sigma(:, :), future_d(:, :)

      if (.not. model%has_priors) then
         info = ssbvar_not_ready
         return
      end if
      horizon = 1
      draws = 1000
      warmup = 500
      rng_seed = 12345
      state_scale = 0.35_dp
      parameter_scale = 0.05_dp
      if (present(h)) horizon = h
      if (present(n_draws)) draws = n_draws
      if (present(burnin)) warmup = burnin
      if (present(seed)) rng_seed = seed
      if (present(sv_state_scale)) state_scale = sv_state_scale
      if (present(sv_parameter_scale)) parameter_scale = sv_parameter_scale
      if (horizon < 1 .or. draws < 1 .or. warmup < 0) then
         info = ssbvar_invalid_input
         return
      end if
      if (state_scale <= 0.0_dp .or. parameter_scale <= 0.0_dp) then
         info = ssbvar_invalid_input
         return
      end if

      call build_future_d(model, horizon, d_pred, future_d, info)
      if (info /= ssbvar_success) return
      model%has_fit = .false.
      call seed_rng(rng_seed)
      model%fit = bvar_fit()

      if (model%priors%sv) then
         call fit_sv_model(model, horizon, draws, warmup, future_d, state_scale, parameter_scale, info)
         if (info /= ssbvar_success) return
         model%has_fit = .true.
         return
      end if

      k = model%setup%k
      p = model%setup%p
      q = model%setup%q
      kp = k * p
      allocate(beta(kp, k), psi(k, q), sigma(k, k))
      beta = model%setup%beta_ols
      psi = model%setup%psi_ols
      sigma = model%setup%sigma_u_ols
      allocate(model%fit%beta(kp, k, draws), model%fit%psi(k, q, draws))
      allocate(model%fit%sigma_u(k, k, draws), model%fit%y_pred(draws, horizon, k))
      allocate(model%fit%d_pred(horizon, q))
      model%fit%d_pred = future_d
      model%fit%n_draws = draws
      model%fit%h = horizon
      model%fit%homoscedastic = .true.
      model%fit%sv_type = ""

      total_iter = warmup + draws
      kept = 0
      do iter = 1, total_iter
         call sample_beta_conditional(model, psi, sigma, beta, info)
         if (info /= ssbvar_success) return
         call sample_psi_conditional(model, beta, sigma, psi, info)
         if (info /= ssbvar_success) return
         call sample_sigma_conditional(model, beta, psi, sigma, info)
         if (info /= ssbvar_success) return
         if (iter > warmup) then
            kept = kept + 1
            model%fit%beta(:, :, kept) = beta
            model%fit%psi(:, :, kept) = psi
            model%fit%sigma_u(:, :, kept) = sigma
            call predictive_draw(model, beta, psi, sigma, future_d, model%fit%y_pred(kept, :, :), info)
            if (info /= ssbvar_success) return
         end if
      end do
      call compute_fit_summaries(model)
      model%has_fit = .true.
      info = ssbvar_success
   end subroutine fit_model


   subroutine fit_sv_model(model, horizon, draws, warmup, future_d, state_scale, parameter_scale, info)
      type(bvar_model), intent(inout) :: model !! BVAR object with validated RW or AR1 stochastic-volatility priors.
      integer, intent(in) :: horizon !! Positive forecast horizon used for generated predictive draws.
      integer, intent(in) :: draws !! Positive number of retained posterior draws.
      integer, intent(in) :: warmup !! Nonnegative number of discarded MCMC iterations.
      real(dp), intent(in) :: future_d(:, :) !! Future deterministic terms with shape (H,q).
      real(dp), intent(in) :: state_scale !! Random-walk proposal standard deviation for latent log-volatility states.
      real(dp), intent(in) :: parameter_scale !! Random-walk proposal standard deviation for AR1 gamma parameters.
      integer, intent(out) :: info !! Status code; zero indicates SV posterior and predictive sampling succeeded.
      integer :: n, k, p, q, kp, total_iter, iter, kept, i, t
      real(dp), allocatable :: beta(:, :), psi(:, :), a(:, :), log_lambda(:, :), phi(:)
      real(dp), allocatable :: gamma_0(:), gamma_1(:), phi_cov(:, :), residuals(:, :), transformed(:, :)
      real(dp), allocatable :: sigma_time(:, :, :)
      real(dp) :: denom

      n = model%setup%n
      k = model%setup%k
      p = model%setup%p
      q = model%setup%q
      kp = k * p
      allocate(beta(kp, k), psi(k, q), a(k, k), log_lambda(n, k))
      allocate(residuals(n, k), transformed(n, k), sigma_time(k, k, n))
      beta = model%setup%beta_ols
      psi = model%setup%psi_ols
      call initialize_a_from_prior(model, a)
      do i = 1, k
         log_lambda(1, i) = model%priors%sv_priors%theta_log_lambda_1(i)
      end do

      if (trim(model%priors%sv_type) == "RW") then
         allocate(phi(k))
         do i = 1, k
            denom = max(model%priors%sv_priors%alpha_phi(i), 1.0e-6_dp)
            phi(i) = model%priors%sv_priors%beta_phi(i) / denom
            phi(i) = max(phi(i), 1.0e-6_dp)
            do t = 2, n
               log_lambda(t, i) = log_lambda(t - 1, i)
            end do
         end do
      else if (trim(model%priors%sv_type) == "AR1") then
         allocate(gamma_0(k), gamma_1(k), phi_cov(k, k))
         gamma_0 = model%priors%sv_priors%theta_gamma_0
         gamma_1 = max(-0.95_dp, min(0.95_dp, model%priors%sv_priors%theta_gamma_1))
         denom = real(max(model%priors%sv_priors%m_phi - k - 1, 1), dp)
         phi_cov = model%priors%sv_priors%v_phi / denom
         do i = 1, k
            phi_cov(i, i) = max(phi_cov(i, i), 1.0e-6_dp)
         end do
         do t = 2, n
            log_lambda(t, :) = gamma_0 + gamma_1 * log_lambda(t - 1, :)
         end do
      else
         info = ssbvar_invalid_input
         return
      end if

      allocate(model%fit%beta(kp, k, draws), model%fit%psi(k, q, draws))
      allocate(model%fit%y_pred(draws, horizon, k), model%fit%d_pred(horizon, q))
      allocate(model%fit%a(k, k, draws), model%fit%log_lambda(n, k, draws))
      allocate(model%fit%sigma_u_time(k, k, n, draws))
      allocate(model%fit%sigma_u_pred(k, k, horizon, draws))
      if (trim(model%priors%sv_type) == "RW") then
         allocate(model%fit%phi(k, draws))
      else
         allocate(model%fit%gamma_0(k, draws), model%fit%gamma_1(k, draws))
         allocate(model%fit%phi_cov(k, k, draws))
      end if
      model%fit%d_pred = future_d
      model%fit%n_draws = draws
      model%fit%h = horizon
      model%fit%homoscedastic = .false.
      model%fit%sv_type = model%priors%sv_type

      total_iter = warmup + draws
      kept = 0
      do iter = 1, total_iter
         call sample_beta_conditional_sv(model, psi, a, log_lambda, beta, info)
         if (info /= ssbvar_success) return
         call sample_psi_conditional_sv(model, beta, a, log_lambda, psi, info)
         if (info /= ssbvar_success) return
         call compute_residuals(model, beta, psi, residuals)
         call sample_a_conditional_sv(model, residuals, log_lambda, a, info)
         if (info /= ssbvar_success) return
         call transform_residuals(a, residuals, transformed)

         if (trim(model%priors%sv_type) == "RW") then
            call update_log_lambda_rw(model, transformed, phi, state_scale, log_lambda)
            call sample_phi_rw(model, log_lambda, phi)
         else
            call update_log_lambda_ar1(model, transformed, gamma_0, gamma_1, phi_cov, &
               state_scale, log_lambda, info)
            if (info /= ssbvar_success) return
            call update_gamma_ar1(model, log_lambda, phi_cov, parameter_scale, gamma_0, gamma_1, info)
            if (info /= ssbvar_success) return
            call sample_phi_cov_ar1(model, log_lambda, gamma_0, gamma_1, phi_cov, info)
            if (info /= ssbvar_success) return
         end if

         if (iter > warmup) then
            kept = kept + 1
            model%fit%beta(:, :, kept) = beta
            model%fit%psi(:, :, kept) = psi
            model%fit%a(:, :, kept) = a
            model%fit%log_lambda(:, :, kept) = log_lambda
            call construct_sigma_time(a, log_lambda, sigma_time, info)
            if (info /= ssbvar_success) return
            model%fit%sigma_u_time(:, :, :, kept) = sigma_time
            if (trim(model%priors%sv_type) == "RW") then
               model%fit%phi(:, kept) = phi
               call predictive_draw_sv(model, beta, psi, a, log_lambda(n, :), future_d, phi=phi, &
                  y_pred=model%fit%y_pred(kept, :, :), sigma_pred=model%fit%sigma_u_pred(:, :, :, kept), info=info)
            else
               model%fit%gamma_0(:, kept) = gamma_0
               model%fit%gamma_1(:, kept) = gamma_1
               model%fit%phi_cov(:, :, kept) = phi_cov
               call predictive_draw_sv(model, beta, psi, a, log_lambda(n, :), future_d, &
                  gamma_0=gamma_0, gamma_1=gamma_1, phi_cov=phi_cov, &
                  y_pred=model%fit%y_pred(kept, :, :), sigma_pred=model%fit%sigma_u_pred(:, :, :, kept), info=info)
            end if
            if (info /= ssbvar_success) return
         end if
      end do
      call compute_fit_summaries(model)
      info = ssbvar_success
   end subroutine fit_sv_model

   pure subroutine initialize_a_from_prior(model, a)
      type(bvar_model), intent(in) :: model !! BVAR object providing the row-major prior means of free A coefficients.
      real(dp), intent(out) :: a(:, :) !! Unit-lower-triangular contemporaneous matrix initialized at prior means.
      integer :: i, j, idx, k

      k = model%setup%k
      a = 0.0_dp
      do i = 1, k
         a(i, i) = 1.0_dp
      end do
      idx = 0
      do i = 2, k
         do j = 1, i - 1
            idx = idx + 1
            a(i, j) = model%priors%sv_priors%theta_a(idx)
         end do
      end do
   end subroutine initialize_a_from_prior

   subroutine sample_beta_conditional_sv(model, psi, a_mat, log_lambda, beta, info)
      type(bvar_model), intent(in) :: model !! Prepared BVAR model supplying data, lag matrices, and beta prior.
      real(dp), intent(in) :: psi(:, :) !! Current steady-state matrix with shape (k,q).
      real(dp), intent(in) :: a_mat(:, :) !! Current unit-lower-triangular contemporaneous matrix A.
      real(dp), intent(in) :: log_lambda(:, :) !! Current latent log volatilities with shape (N,k).
      real(dp), intent(out) :: beta(:, :) !! Newly sampled autoregressive coefficient matrix with shape (k*p,k).
      integer, intent(out) :: info !! Status code; zero indicates the time-varying Gaussian conditional draw succeeded.
      integer :: n, k, p, q, kp, d_beta, t, lag, i, j, ii, jj, idx1, idx2
      real(dp), allocatable :: bmat(:, :), y_center(:, :), precision(:, :), rhs(:), mean_vec(:)
      real(dp), allocatable :: covariance(:, :), draw(:), basis(:, :), sigma_inv(:, :), weighted_y(:)
      real(dp) :: logdet, value

      n = model%setup%n
      k = model%setup%k
      p = model%setup%p
      q = model%setup%q
      kp = k * p
      d_beta = kp * k
      allocate(bmat(n, kp), y_center(n, k), basis(q * p, kp))
      basis = 0.0_dp
      do lag = 1, p
         basis((lag - 1) * q + 1:lag * q, (lag - 1) * k + 1:lag * k) = transpose(psi)
      end do
      bmat = model%setup%w - matmul(model%setup%q_lag, basis)
      y_center = model%setup%y - matmul(model%setup%d, transpose(psi))
      allocate(precision(d_beta, d_beta), rhs(d_beta), mean_vec(d_beta))
      allocate(sigma_inv(k, k), weighted_y(k))
      precision = 0.0_dp
      rhs = 0.0_dp
      do i = 1, d_beta
         precision(i, i) = 1.0_dp / model%priors%omega_beta(i, i)
         rhs(i) = model%priors%theta_beta(i) / model%priors%omega_beta(i, i)
      end do
      do t = 1, n
         call sv_sigma_inverse(a_mat, log_lambda(t, :), sigma_inv)
         weighted_y = matmul(sigma_inv, y_center(t, :))
         do i = 1, k
            do ii = 1, kp
               idx1 = (i - 1) * kp + ii
               rhs(idx1) = rhs(idx1) + bmat(t, ii) * weighted_y(i)
               do j = 1, k
                  do jj = 1, kp
                     idx2 = (j - 1) * kp + jj
                     value = sigma_inv(i, j) * bmat(t, ii) * bmat(t, jj)
                     precision(idx1, idx2) = precision(idx1, idx2) + value
                  end do
               end do
            end do
         end do
      end do
      call solve_spd(precision, rhs, mean_vec, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      call spd_inverse_logdet(precision, covariance, logdet, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      allocate(draw(d_beta))
      call random_multivariate_normal(mean_vec, covariance, draw, info)
      if (info /= ssbvar_success) return
      beta = reshape(draw, [kp, k])
      info = ssbvar_success
   end subroutine sample_beta_conditional_sv

   subroutine sample_psi_conditional_sv(model, beta, a_mat, log_lambda, psi, info)
      type(bvar_model), intent(in) :: model !! Prepared BVAR model supplying deterministic terms and the Psi prior.
      real(dp), intent(in) :: beta(:, :) !! Current autoregressive coefficient matrix with shape (k*p,k).
      real(dp), intent(in) :: a_mat(:, :) !! Current unit-lower-triangular contemporaneous matrix A.
      real(dp), intent(in) :: log_lambda(:, :) !! Current latent log volatilities with shape (N,k).
      real(dp), intent(out) :: psi(:, :) !! Newly sampled steady-state matrix with shape (k,q).
      integer, intent(out) :: info !! Status code; zero indicates the time-varying Gaussian conditional draw succeeded.
      integer :: n, k, p, q, kp, d_psi, t, lag, i, b, idx
      real(dp), allocatable :: precision(:, :), rhs(:), mean_vec(:), covariance(:, :)
      real(dp), allocatable :: hmat(:, :), y0(:), pi_l(:, :, :), draw(:), temp(:), sigma_inv(:, :)
      real(dp) :: logdet

      n = model%setup%n
      k = model%setup%k
      p = model%setup%p
      q = model%setup%q
      kp = k * p
      d_psi = k * q
      allocate(precision(d_psi, d_psi), rhs(d_psi), mean_vec(d_psi))
      allocate(hmat(k, d_psi), y0(k), pi_l(k, k, p), temp(k), sigma_inv(k, k))
      precision = 0.0_dp
      rhs = 0.0_dp
      do idx = 1, d_psi
         precision(idx, idx) = 1.0_dp / model%priors%omega_psi(idx, idx)
         rhs(idx) = model%priors%theta_psi(idx) / model%priors%omega_psi(idx, idx)
      end do
      do lag = 1, p
         pi_l(:, :, lag) = transpose(beta((lag - 1) * k + 1:lag * k, :))
      end do
      do t = 1, n
         hmat = 0.0_dp
         do b = 1, q
            do i = 1, k
               idx = (b - 1) * k + i
               hmat(i, idx) = model%setup%d(t, b)
               do lag = 1, p
                  hmat(:, idx) = hmat(:, idx) - &
                     model%setup%q_lag(t, (lag - 1) * q + b) * pi_l(:, i, lag)
               end do
            end do
         end do
         y0 = model%setup%y(t, :) - matmul(model%setup%w(t, :), beta)
         call sv_sigma_inverse(a_mat, log_lambda(t, :), sigma_inv)
         precision = precision + matmul(transpose(hmat), matmul(sigma_inv, hmat))
         temp = matmul(sigma_inv, y0)
         rhs = rhs + matmul(transpose(hmat), temp)
      end do
      call solve_spd(precision, rhs, mean_vec, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      call spd_inverse_logdet(precision, covariance, logdet, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      allocate(draw(d_psi))
      call random_multivariate_normal(mean_vec, covariance, draw, info)
      if (info /= ssbvar_success) return
      psi = reshape(draw, [k, q])
      info = ssbvar_success
   end subroutine sample_psi_conditional_sv

   pure subroutine compute_residuals(model, beta, psi, residuals)
      type(bvar_model), intent(in) :: model !! Prepared BVAR model containing observed and lagged regressors.
      real(dp), intent(in) :: beta(:, :) !! Autoregressive coefficient matrix with shape (k*p,k).
      real(dp), intent(in) :: psi(:, :) !! Steady-state coefficient matrix with shape (k,q).
      real(dp), intent(out) :: residuals(:, :) !! Reduced-form residual matrix with shape (N,k).
      integer :: p, q, k, lag
      real(dp), allocatable :: basis(:, :), bmat(:, :)

      p = model%setup%p
      q = model%setup%q
      k = model%setup%k
      allocate(basis(q * p, k * p), bmat(model%setup%n, k * p))
      basis = 0.0_dp
      do lag = 1, p
         basis((lag - 1) * q + 1:lag * q, (lag - 1) * k + 1:lag * k) = transpose(psi)
      end do
      bmat = model%setup%w - matmul(model%setup%q_lag, basis)
      residuals = model%setup%y - matmul(model%setup%d, transpose(psi)) - matmul(bmat, beta)
   end subroutine compute_residuals

   subroutine sample_a_conditional_sv(model, residuals, log_lambda, a_mat, info)
      type(bvar_model), intent(in) :: model !! BVAR object providing the independent normal priors for free A coefficients.
      real(dp), intent(in) :: residuals(:, :) !! Current reduced-form residuals with shape (N,k).
      real(dp), intent(in) :: log_lambda(:, :) !! Current latent log volatilities with shape (N,k).
      real(dp), intent(inout) :: a_mat(:, :) !! Unit-lower-triangular A matrix updated row by row.
      integer, intent(out) :: info !! Status code; zero indicates all free A rows were sampled successfully.
      integer :: n, k, row, d, start_idx, i, j, t
      real(dp), allocatable :: precision(:, :), rhs(:), mean_vec(:), covariance(:, :), draw(:), x(:)
      real(dp) :: weight, y_value, logdet

      n = model%setup%n
      k = model%setup%k
      do i = 1, k
         a_mat(i, i) = 1.0_dp
         if (i < k) a_mat(i, i + 1:k) = 0.0_dp
      end do
      do row = 2, k
         d = row - 1
         start_idx = 1 + (row - 2) * (row - 1) / 2
         allocate(precision(d, d), rhs(d), mean_vec(d), x(d))
         precision = 0.0_dp
         rhs = 0.0_dp
         do i = 1, d
            precision(i, i) = 1.0_dp / model%priors%sv_priors%omega_a(start_idx + i - 1, start_idx + i - 1)
            rhs(i) = model%priors%sv_priors%theta_a(start_idx + i - 1) / &
               model%priors%sv_priors%omega_a(start_idx + i - 1, start_idx + i - 1)
         end do
         do t = 1, n
            x = residuals(t, 1:d)
            y_value = -residuals(t, row)
            weight = exp(-log_lambda(t, row))
            do i = 1, d
               rhs(i) = rhs(i) + weight * x(i) * y_value
               do j = 1, d
                  precision(i, j) = precision(i, j) + weight * x(i) * x(j)
               end do
            end do
         end do
         call solve_spd(precision, rhs, mean_vec, info)
         if (info /= 0) then
            info = ssbvar_linalg_failure
            return
         end if
         call spd_inverse_logdet(precision, covariance, logdet, info)
         if (info /= 0) then
            info = ssbvar_linalg_failure
            return
         end if
         allocate(draw(d))
         call random_multivariate_normal(mean_vec, covariance, draw, info)
         if (info /= ssbvar_success) return
         a_mat(row, 1:d) = draw
         deallocate(precision, rhs, mean_vec, x, covariance, draw)
      end do
      info = ssbvar_success
   end subroutine sample_a_conditional_sv

   pure subroutine transform_residuals(a_mat, residuals, transformed)
      real(dp), intent(in) :: a_mat(:, :) !! Unit-lower-triangular contemporaneous matrix A.
      real(dp), intent(in) :: residuals(:, :) !! Reduced-form residual matrix with rows as time.
      real(dp), intent(out) :: transformed(:, :) !! Structural residuals A*u_t with the same time-by-variable shape.
      integer :: t

      do t = 1, size(residuals, 1)
         transformed(t, :) = matmul(a_mat, residuals(t, :))
      end do
   end subroutine transform_residuals

   pure subroutine sv_sigma_inverse(a_mat, log_lambda_t, sigma_inv)
      real(dp), intent(in) :: a_mat(:, :) !! Unit-lower-triangular contemporaneous matrix A.
      real(dp), intent(in) :: log_lambda_t(:) !! Log structural variances for one time point, length k.
      real(dp), intent(out) :: sigma_inv(:, :) !! Inverse reduced-form covariance A' diag(exp(-log_lambda)) A.
      real(dp) :: scaled(size(a_mat, 1), size(a_mat, 2))
      integer :: i

      scaled = a_mat
      do i = 1, size(a_mat, 1)
         scaled(i, :) = exp(-log_lambda_t(i)) * a_mat(i, :)
      end do
      sigma_inv = matmul(transpose(a_mat), scaled)
   end subroutine sv_sigma_inverse

   subroutine update_log_lambda_rw(model, transformed, phi, proposal_scale, log_lambda)
      type(bvar_model), intent(in) :: model !! BVAR object providing RW initial-state priors.
      real(dp), intent(in) :: transformed(:, :) !! Structural residuals with shape (N,k).
      real(dp), intent(in) :: phi(:) !! Positive RW innovation variances for each log-volatility series.
      real(dp), intent(in) :: proposal_scale !! Positive scalar random-walk proposal standard deviation.
      real(dp), intent(inout) :: log_lambda(:, :) !! Latent log-volatility paths updated by scalar Metropolis steps.
      integer :: i, t
      real(dp) :: proposal, old_lp, new_lp, u
      real(dp) :: z(1)

      do i = 1, size(log_lambda, 2)
         do t = 1, size(log_lambda, 1)
            call random_standard_normal_vector(z)
            proposal = log_lambda(t, i) + proposal_scale * z(1)
            old_lp = rw_state_log_kernel(model, transformed, phi, log_lambda, i, t, log_lambda(t, i))
            new_lp = rw_state_log_kernel(model, transformed, phi, log_lambda, i, t, proposal)
            call random_number(u)
            if (log(max(u, tiny(1.0_dp))) < new_lp - old_lp) log_lambda(t, i) = proposal
         end do
      end do
   end subroutine update_log_lambda_rw

   pure real(dp) function rw_state_log_kernel(model, transformed, phi, log_lambda, i, t, candidate) result(value)
      type(bvar_model), intent(in) :: model !! BVAR object supplying the RW initial-state normal prior.
      real(dp), intent(in) :: transformed(:, :) !! Structural residuals with shape (N,k).
      real(dp), intent(in) :: phi(:) !! Positive RW innovation variances, length k.
      real(dp), intent(in) :: log_lambda(:, :) !! Current latent path used for neighboring states.
      integer, intent(in) :: i !! One-based stochastic-volatility series index.
      integer, intent(in) :: t !! One-based time index of the proposed latent state.
      real(dp), intent(in) :: candidate !! Candidate log variance at series i and time t.
      integer :: n
      real(dp) :: diff, variance

      n = size(log_lambda, 1)
      value = -0.5_dp * (candidate + transformed(t, i) ** 2 * exp(-candidate))
      if (t == 1) then
         variance = model%priors%sv_priors%omega_log_lambda_1(i, i)
         diff = candidate - model%priors%sv_priors%theta_log_lambda_1(i)
         value = value - 0.5_dp * diff * diff / variance
      else
         diff = candidate - log_lambda(t - 1, i)
         value = value - 0.5_dp * diff * diff / phi(i)
      end if
      if (t < n) then
         diff = log_lambda(t + 1, i) - candidate
         value = value - 0.5_dp * diff * diff / phi(i)
      end if
   end function rw_state_log_kernel

   subroutine sample_phi_rw(model, log_lambda, phi)
      type(bvar_model), intent(in) :: model !! BVAR object supplying inverse-gamma priors for RW innovation variances.
      real(dp), intent(in) :: log_lambda(:, :) !! Current latent RW log-volatility paths with shape (N,k).
      real(dp), intent(out) :: phi(:) !! Newly sampled positive RW innovation variances, length k.
      integer :: i, n
      real(dp) :: shape, scale, gamma_value

      n = size(log_lambda, 1)
      do i = 1, size(log_lambda, 2)
         shape = model%priors%sv_priors%alpha_phi(i) + 0.5_dp * real(n - 1, dp)
         scale = model%priors%sv_priors%beta_phi(i)
         if (n > 1) scale = scale + 0.5_dp * sum((log_lambda(2:n, i) - log_lambda(1:n - 1, i)) ** 2)
         call random_gamma(shape, gamma_value)
         phi(i) = scale / gamma_value
      end do
   end subroutine sample_phi_rw

   subroutine update_log_lambda_ar1(model, transformed, gamma_0, gamma_1, phi_cov, proposal_scale, &
      log_lambda, info)
      type(bvar_model), intent(in) :: model !! BVAR object providing the AR1 initial-state prior.
      real(dp), intent(in) :: transformed(:, :) !! Structural residuals with shape (N,k).
      real(dp), intent(in) :: gamma_0(:) !! Current AR1 log-volatility intercepts, length k.
      real(dp), intent(in) :: gamma_1(:) !! Current AR1 log-volatility slopes constrained to (-1,1), length k.
      real(dp), intent(in) :: phi_cov(:, :) !! Current covariance of AR1 log-volatility innovations.
      real(dp), intent(in) :: proposal_scale !! Positive scalar random-walk proposal standard deviation.
      real(dp), intent(inout) :: log_lambda(:, :) !! Latent AR1 log-volatility paths updated by scalar Metropolis steps.
      integer, intent(out) :: info !! Status code; zero indicates covariance inversion and all state updates succeeded.
      integer :: i, t
      real(dp) :: proposal, old_lp, new_lp, u, logdet
      real(dp), allocatable :: phi_inv(:, :)
      real(dp) :: z(1)

      call spd_inverse_logdet(phi_cov, phi_inv, logdet, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      do i = 1, size(log_lambda, 2)
         do t = 1, size(log_lambda, 1)
            call random_standard_normal_vector(z)
            proposal = log_lambda(t, i) + proposal_scale * z(1)
            old_lp = ar1_state_log_kernel(model, transformed, gamma_0, gamma_1, phi_inv, &
               log_lambda, i, t, log_lambda(t, i))
            new_lp = ar1_state_log_kernel(model, transformed, gamma_0, gamma_1, phi_inv, &
               log_lambda, i, t, proposal)
            call random_number(u)
            if (log(max(u, tiny(1.0_dp))) < new_lp - old_lp) log_lambda(t, i) = proposal
         end do
      end do
      info = ssbvar_success
   end subroutine update_log_lambda_ar1

   pure real(dp) function ar1_state_log_kernel(model, transformed, gamma_0, gamma_1, phi_inv, &
      log_lambda, i, t, candidate) result(value)
      type(bvar_model), intent(in) :: model !! BVAR object supplying the initial log-volatility prior.
      real(dp), intent(in) :: transformed(:, :) !! Structural residuals with shape (N,k).
      real(dp), intent(in) :: gamma_0(:) !! Current AR1 intercept vector.
      real(dp), intent(in) :: gamma_1(:) !! Current AR1 slope vector.
      real(dp), intent(in) :: phi_inv(:, :) !! Inverse covariance of AR1 log-volatility innovations.
      real(dp), intent(in) :: log_lambda(:, :) !! Current latent path used for neighboring states.
      integer, intent(in) :: i !! One-based stochastic-volatility series index.
      integer, intent(in) :: t !! One-based time index of the proposed latent state.
      real(dp), intent(in) :: candidate !! Candidate log variance at series i and time t.
      integer :: n
      real(dp) :: diff, variance
      real(dp) :: residual(size(log_lambda, 2))

      n = size(log_lambda, 1)
      value = -0.5_dp * (candidate + transformed(t, i) ** 2 * exp(-candidate))
      if (t == 1) then
         variance = model%priors%sv_priors%omega_log_lambda_1(i, i)
         diff = candidate - model%priors%sv_priors%theta_log_lambda_1(i)
         value = value - 0.5_dp * diff * diff / variance
      else
         residual = log_lambda(t, :) - gamma_0 - gamma_1 * log_lambda(t - 1, :)
         residual(i) = candidate - gamma_0(i) - gamma_1(i) * log_lambda(t - 1, i)
         value = value - 0.5_dp * dot_product(residual, matmul(phi_inv, residual))
      end if
      if (t < n) then
         residual = log_lambda(t + 1, :) - gamma_0 - gamma_1 * log_lambda(t, :)
         residual(i) = log_lambda(t + 1, i) - gamma_0(i) - gamma_1(i) * candidate
         value = value - 0.5_dp * dot_product(residual, matmul(phi_inv, residual))
      end if
   end function ar1_state_log_kernel

   subroutine update_gamma_ar1(model, log_lambda, phi_cov, proposal_scale, gamma_0, gamma_1, info)
      type(bvar_model), intent(in) :: model !! BVAR object supplying independent normal priors for AR1 gamma parameters.
      real(dp), intent(in) :: log_lambda(:, :) !! Current latent AR1 log-volatility paths with shape (N,k).
      real(dp), intent(in) :: phi_cov(:, :) !! Current covariance of AR1 log-volatility innovations.
      real(dp), intent(in) :: proposal_scale !! Positive random-walk proposal standard deviation for gamma parameters.
      real(dp), intent(inout) :: gamma_0(:) !! AR1 intercept vector updated by scalar Metropolis steps.
      real(dp), intent(inout) :: gamma_1(:) !! AR1 slope vector updated within the open stationary interval (-1,1).
      integer, intent(out) :: info !! Status code; zero indicates covariance inversion and all gamma updates succeeded.
      integer :: i
      real(dp) :: proposal, old_lp, new_lp, u, logdet
      real(dp), allocatable :: phi_inv(:, :), work_0(:), work_1(:)
      real(dp) :: z(1)

      call spd_inverse_logdet(phi_cov, phi_inv, logdet, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      allocate(work_0(size(gamma_0)), work_1(size(gamma_1)))
      do i = 1, size(gamma_0)
         call random_standard_normal_vector(z)
         proposal = gamma_0(i) + proposal_scale * z(1)
         old_lp = ar1_transition_log_kernel(log_lambda, gamma_0, gamma_1, phi_inv) - &
            0.5_dp * (gamma_0(i) - model%priors%sv_priors%theta_gamma_0(i)) ** 2 / &
            model%priors%sv_priors%omega_gamma_0(i, i)
         work_0 = gamma_0
         work_0(i) = proposal
         new_lp = ar1_transition_log_kernel(log_lambda, work_0, gamma_1, phi_inv) - &
            0.5_dp * (proposal - model%priors%sv_priors%theta_gamma_0(i)) ** 2 / &
            model%priors%sv_priors%omega_gamma_0(i, i)
         call random_number(u)
         if (log(max(u, tiny(1.0_dp))) < new_lp - old_lp) gamma_0(i) = proposal

         call random_standard_normal_vector(z)
         proposal = gamma_1(i) + proposal_scale * z(1)
         if (abs(proposal) < 1.0_dp) then
            old_lp = ar1_transition_log_kernel(log_lambda, gamma_0, gamma_1, phi_inv) - &
               0.5_dp * (gamma_1(i) - model%priors%sv_priors%theta_gamma_1(i)) ** 2 / &
               model%priors%sv_priors%omega_gamma_1(i, i)
            work_1 = gamma_1
            work_1(i) = proposal
            new_lp = ar1_transition_log_kernel(log_lambda, gamma_0, work_1, phi_inv) - &
               0.5_dp * (proposal - model%priors%sv_priors%theta_gamma_1(i)) ** 2 / &
               model%priors%sv_priors%omega_gamma_1(i, i)
            call random_number(u)
            if (log(max(u, tiny(1.0_dp))) < new_lp - old_lp) gamma_1(i) = proposal
         end if
      end do
      info = ssbvar_success
   end subroutine update_gamma_ar1

   pure real(dp) function ar1_transition_log_kernel(log_lambda, gamma_0, gamma_1, phi_inv) result(value)
      real(dp), intent(in) :: log_lambda(:, :) !! Latent AR1 log-volatility paths with shape (N,k).
      real(dp), intent(in) :: gamma_0(:) !! Candidate AR1 intercept vector.
      real(dp), intent(in) :: gamma_1(:) !! Candidate AR1 slope vector.
      real(dp), intent(in) :: phi_inv(:, :) !! Inverse covariance of AR1 log-volatility innovations.
      integer :: t
      real(dp) :: residual(size(log_lambda, 2))

      value = 0.0_dp
      do t = 2, size(log_lambda, 1)
         residual = log_lambda(t, :) - gamma_0 - gamma_1 * log_lambda(t - 1, :)
         value = value - 0.5_dp * dot_product(residual, matmul(phi_inv, residual))
      end do
   end function ar1_transition_log_kernel

   subroutine sample_phi_cov_ar1(model, log_lambda, gamma_0, gamma_1, phi_cov, info)
      type(bvar_model), intent(in) :: model !! BVAR object supplying the inverse-Wishart prior for AR1 innovation covariance Phi.
      real(dp), intent(in) :: log_lambda(:, :) !! Current latent AR1 log-volatility paths with shape (N,k).
      real(dp), intent(in) :: gamma_0(:) !! Current AR1 intercept vector, length k.
      real(dp), intent(in) :: gamma_1(:) !! Current AR1 slope vector, length k.
      real(dp), intent(out) :: phi_cov(:, :) !! Newly sampled AR1 innovation covariance matrix Phi.
      integer, intent(out) :: info !! Status code; zero indicates the inverse-Wishart draw succeeded.
      integer :: t, n, k, degrees
      real(dp), allocatable :: scale(:, :)
      real(dp) :: residual(size(log_lambda, 2))

      n = size(log_lambda, 1)
      k = size(log_lambda, 2)
      allocate(scale(k, k))
      scale = model%priors%sv_priors%v_phi
      do t = 2, n
         residual = log_lambda(t, :) - gamma_0 - gamma_1 * log_lambda(t - 1, :)
         scale = scale + spread(residual, 2, k) * spread(residual, 1, k)
      end do
      degrees = model%priors%sv_priors%m_phi + n - 1
      call random_inverse_wishart(degrees, scale, phi_cov, info)
   end subroutine sample_phi_cov_ar1

   subroutine construct_sigma_time(a_mat, log_lambda, sigma_time, info)
      real(dp), intent(in) :: a_mat(:, :) !! Unit-lower-triangular contemporaneous matrix A.
      real(dp), intent(in) :: log_lambda(:, :) !! Latent log structural variances with shape (N,k).
      real(dp), intent(out) :: sigma_time(:, :, :) !! Reduced-form covariance sequence with shape (k,k,N).
      integer, intent(out) :: info !! Status code; zero indicates inversion and all covariance constructions succeeded.
      integer :: t, j, k
      real(dp), allocatable :: a_inv(:, :), scaled(:, :)

      k = size(a_mat, 1)
      call inverse_matrix(a_mat, a_inv, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      allocate(scaled(k, k))
      do t = 1, size(log_lambda, 1)
         scaled = a_inv
         do j = 1, k
            scaled(:, j) = exp(log_lambda(t, j)) * a_inv(:, j)
         end do
         sigma_time(:, :, t) = matmul(scaled, transpose(a_inv))
      end do
      info = ssbvar_success
   end subroutine construct_sigma_time

   subroutine predictive_draw_sv(model, beta, psi, a_mat, last_log_lambda, d_pred, phi, gamma_0, gamma_1, &
      phi_cov, y_pred, sigma_pred, info)
      type(bvar_model), intent(in) :: model !! Fitted-model setup supplying observed lags for predictive recursion.
      real(dp), intent(in) :: beta(:, :) !! Autoregressive coefficient draw with shape (k*p,k).
      real(dp), intent(in) :: psi(:, :) !! Steady-state parameter draw with shape (k,q).
      real(dp), intent(in) :: a_mat(:, :) !! Unit-lower-triangular contemporaneous matrix A.
      real(dp), intent(in) :: last_log_lambda(:) !! Last in-sample latent log-volatility state, length k.
      real(dp), intent(in) :: d_pred(:, :) !! Future deterministic terms with shape (H,q).
      real(dp), intent(in), optional :: phi(:) !! RW log-volatility innovation variances, required for RW models.
      real(dp), intent(in), optional :: gamma_0(:) !! AR1 log-volatility intercepts, required for AR1 models.
      real(dp), intent(in), optional :: gamma_1(:) !! AR1 log-volatility slopes, required for AR1 models.
      real(dp), intent(in), optional :: phi_cov(:, :) !! AR1 innovation covariance Phi, required for AR1 models.
      real(dp), intent(out) :: y_pred(:, :) !! Joint predictive draw with shape (H,k).
      real(dp), intent(out) :: sigma_pred(:, :, :) !! Future reduced-form covariances with shape (k,k,H).
      integer, intent(out) :: info !! Status code; zero indicates all stochastic-volatility forecast draws were generated.
      integer :: h, horizon, k, p, lag, n, j
      real(dp), allocatable :: pi_l(:, :, :), yhat(:), innovation(:), mu_hist(:), mu_fore(:)
      real(dp), allocatable :: a_inv(:, :), l_factor(:, :), current_log(:), z(:), nu(:), zero(:)

      horizon = size(d_pred, 1)
      k = model%setup%k
      p = model%setup%p
      n = model%setup%n
      allocate(pi_l(k, k, p), yhat(k), innovation(k), mu_hist(k), mu_fore(k))
      allocate(current_log(k), z(k), nu(k), zero(k))
      zero = 0.0_dp
      current_log = last_log_lambda
      call inverse_matrix(a_mat, a_inv, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      do lag = 1, p
         pi_l(:, :, lag) = transpose(beta((lag - 1) * k + 1:lag * k, :))
      end do
      do h = 1, horizon
         if (trim(model%priors%sv_type) == "RW") then
            if (.not. present(phi)) then
               info = ssbvar_invalid_input
               return
            end if
            call random_standard_normal_vector(z)
            current_log = current_log + sqrt(phi) * z
         else
            if (.not. present(gamma_0) .or. .not. present(gamma_1) .or. .not. present(phi_cov)) then
               info = ssbvar_invalid_input
               return
            end if
            call random_multivariate_normal(zero, phi_cov, nu, info)
            if (info /= ssbvar_success) return
            current_log = gamma_0 + gamma_1 * current_log + nu
         end if
         allocate(l_factor(k, k))
         l_factor = a_inv
         do j = 1, k
            l_factor(:, j) = exp(0.5_dp * current_log(j)) * a_inv(:, j)
         end do
         sigma_pred(:, :, h) = matmul(l_factor, transpose(l_factor))
         call random_standard_normal_vector(z)
         innovation = matmul(l_factor, z)
         deallocate(l_factor)

         yhat = matmul(psi, d_pred(h, :))
         do lag = 1, min(h - 1, p)
            mu_fore = matmul(psi, d_pred(h - lag, :))
            yhat = yhat + matmul(pi_l(:, :, lag), y_pred(h - lag, :) - mu_fore)
         end do
         if (h <= p) then
            do lag = h, p
               mu_hist = matmul(psi, model%setup%x(n + h - lag, :))
               yhat = yhat + matmul(pi_l(:, :, lag), model%setup%y(n + h - lag, :) - mu_hist)
            end do
         end if
         y_pred(h, :) = yhat + innovation
      end do
      info = ssbvar_success
   end subroutine predictive_draw_sv

   subroutine sample_beta_conditional(model, psi, sigma, beta, info)
      type(bvar_model), intent(in) :: model !! Prepared BVAR model supplying data, lag matrices, and beta prior.
      real(dp), intent(in) :: psi(:, :) !! Current steady-state matrix with shape (k,q).
      real(dp), intent(in) :: sigma(:, :) !! Current innovation covariance with shape (k,k).
      real(dp), intent(out) :: beta(:, :) !! Newly sampled autoregressive coefficient matrix with shape (k*p,k).
      integer, intent(out) :: info !! Status code; zero indicates the Gaussian conditional draw succeeded.
      integer :: n, k, p, q, kp, d_beta, lag, i, j, a, b, idx1, idx2
      real(dp), allocatable :: bmat(:, :), y_center(:, :), sigma_inv(:, :), xtx(:, :)
      real(dp), allocatable :: precision(:, :), rhs(:), mean_vec(:), covariance(:, :), draw(:)
      real(dp), allocatable :: temp(:, :), basis(:, :)
      real(dp) :: logdet

      n = model%setup%n
      k = model%setup%k
      p = model%setup%p
      q = model%setup%q
      kp = k * p
      d_beta = kp * k
      allocate(bmat(n, kp), y_center(n, k), basis(q * p, kp))
      basis = 0.0_dp
      do lag = 1, p
         basis((lag - 1) * q + 1:lag * q, (lag - 1) * k + 1:lag * k) = transpose(psi)
      end do
      bmat = model%setup%w - matmul(model%setup%q_lag, basis)
      y_center = model%setup%y - matmul(model%setup%d, transpose(psi))

      call spd_inverse_logdet(sigma, sigma_inv, logdet, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      allocate(xtx(kp, kp), precision(d_beta, d_beta), rhs(d_beta), mean_vec(d_beta))
      allocate(temp(kp, k))
      xtx = matmul(transpose(bmat), bmat)
      precision = 0.0_dp
      rhs = 0.0_dp
      do i = 1, d_beta
         precision(i, i) = 1.0_dp / model%priors%omega_beta(i, i)
         rhs(i) = model%priors%theta_beta(i) / model%priors%omega_beta(i, i)
      end do
      do j = 1, k
         do i = 1, k
            do b = 1, kp
               idx2 = (j - 1) * kp + b
               do a = 1, kp
                  idx1 = (i - 1) * kp + a
                  precision(idx1, idx2) = precision(idx1, idx2) + sigma_inv(i, j) * xtx(a, b)
               end do
            end do
         end do
      end do
      temp = matmul(transpose(bmat), matmul(y_center, sigma_inv))
      rhs = rhs + reshape(temp, [d_beta])
      call solve_spd(precision, rhs, mean_vec, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      call spd_inverse_logdet(precision, covariance, logdet, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      allocate(draw(d_beta))
      call random_multivariate_normal(mean_vec, covariance, draw, info)
      if (info /= ssbvar_success) return
      beta = reshape(draw, [kp, k])
      info = ssbvar_success
   end subroutine sample_beta_conditional

   subroutine sample_psi_conditional(model, beta, sigma, psi, info)
      type(bvar_model), intent(in) :: model !! Prepared BVAR model supplying deterministic terms and Psi prior.
      real(dp), intent(in) :: beta(:, :) !! Current autoregressive coefficient matrix with shape (k*p,k).
      real(dp), intent(in) :: sigma(:, :) !! Current innovation covariance matrix with shape (k,k).
      real(dp), intent(out) :: psi(:, :) !! Newly sampled steady-state matrix with shape (k,q).
      integer, intent(out) :: info !! Status code; zero indicates the Gaussian conditional draw succeeded.
      integer :: n, k, p, q, kp, d_psi, t, lag, a, b, idx
      real(dp), allocatable :: sigma_inv(:, :), precision(:, :), rhs(:), mean_vec(:), covariance(:, :)
      real(dp), allocatable :: hmat(:, :), y0(:), pi_l(:, :, :), draw(:), temp(:)
      real(dp) :: logdet

      n = model%setup%n
      k = model%setup%k
      p = model%setup%p
      q = model%setup%q
      kp = k * p
      d_psi = k * q
      call spd_inverse_logdet(sigma, sigma_inv, logdet, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      allocate(precision(d_psi, d_psi), rhs(d_psi), mean_vec(d_psi))
      allocate(hmat(k, d_psi), y0(k), pi_l(k, k, p), temp(k))
      precision = 0.0_dp
      rhs = 0.0_dp
      do idx = 1, d_psi
         precision(idx, idx) = 1.0_dp / model%priors%omega_psi(idx, idx)
         rhs(idx) = model%priors%theta_psi(idx) / model%priors%omega_psi(idx, idx)
      end do
      do lag = 1, p
         pi_l(:, :, lag) = transpose(beta((lag - 1) * k + 1:lag * k, :))
      end do

      do t = 1, n
         hmat = 0.0_dp
         do b = 1, q
            do a = 1, k
               idx = (b - 1) * k + a
               hmat(a, idx) = model%setup%d(t, b)
               do lag = 1, p
                  hmat(:, idx) = hmat(:, idx) - &
                     model%setup%q_lag(t, (lag - 1) * q + b) * pi_l(:, a, lag)
               end do
            end do
         end do
         y0 = model%setup%y(t, :) - matmul(model%setup%w(t, :), beta)
         precision = precision + matmul(transpose(hmat), matmul(sigma_inv, hmat))
         temp = matmul(sigma_inv, y0)
         rhs = rhs + matmul(transpose(hmat), temp)
      end do
      call solve_spd(precision, rhs, mean_vec, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      call spd_inverse_logdet(precision, covariance, logdet, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      allocate(draw(d_psi))
      call random_multivariate_normal(mean_vec, covariance, draw, info)
      if (info /= ssbvar_success) return
      psi = reshape(draw, [k, q])
      info = ssbvar_success
   end subroutine sample_psi_conditional

   subroutine sample_sigma_conditional(model, beta, psi, sigma, info)
      type(bvar_model), intent(in) :: model !! Prepared BVAR model supplying observations and covariance prior choice.
      real(dp), intent(in) :: beta(:, :) !! Current autoregressive coefficient matrix with shape (k*p,k).
      real(dp), intent(in) :: psi(:, :) !! Current steady-state matrix with shape (k,q).
      real(dp), intent(out) :: sigma(:, :) !! Newly sampled innovation covariance matrix with shape (k,k).
      integer, intent(out) :: info !! Status code; zero indicates the inverse-Wishart draw succeeded.
      integer :: n, k, p, q, lag
      real(dp), allocatable :: basis(:, :), bmat(:, :), residuals(:, :), scale(:, :)
      integer :: degrees

      n = model%setup%n
      k = model%setup%k
      p = model%setup%p
      q = model%setup%q
      allocate(basis(q * p, k * p), bmat(n, k * p), residuals(n, k), scale(k, k))
      basis = 0.0_dp
      do lag = 1, p
         basis((lag - 1) * q + 1:lag * q, (lag - 1) * k + 1:lag * k) = transpose(psi)
      end do
      bmat = model%setup%w - matmul(model%setup%q_lag, basis)
      residuals = model%setup%y - matmul(model%setup%d, transpose(psi)) - matmul(bmat, beta)
      scale = matmul(transpose(residuals), residuals)
      if (model%priors%jeffreys) then
         degrees = n
      else
         degrees = model%priors%m + n
         scale = scale + model%priors%v
      end if
      if (degrees <= k - 1) then
         info = ssbvar_invalid_input
         return
      end if
      call random_inverse_wishart(degrees, scale, sigma, info)
   end subroutine sample_sigma_conditional

   subroutine predictive_draw(model, beta, psi, sigma, d_pred, y_pred, info)
      type(bvar_model), intent(in) :: model !! Fitted model and observed sample used to initialize forecast lags.
      real(dp), intent(in) :: beta(:, :) !! Autoregressive coefficient draw with shape (k*p,k).
      real(dp), intent(in) :: psi(:, :) !! Steady-state parameter draw with shape (k,q).
      real(dp), intent(in) :: sigma(:, :) !! Homoscedastic covariance draw with shape (k,k).
      real(dp), intent(in) :: d_pred(:, :) !! Future deterministic values with shape (H,q).
      real(dp), intent(out) :: y_pred(:, :) !! One joint predictive draw with shape (H,k).
      integer, intent(out) :: info !! Status code; zero indicates all forecast innovations were generated.
      integer :: h, horizon, k, p, lag, n
      real(dp), allocatable :: pi_l(:, :, :), yhat(:), innovation(:), mu_hist(:), mu_fore(:)

      horizon = size(d_pred, 1)
      k = model%setup%k
      p = model%setup%p
      n = model%setup%n
      allocate(pi_l(k, k, p), yhat(k), innovation(k), mu_hist(k), mu_fore(k))
      do lag = 1, p
         pi_l(:, :, lag) = transpose(beta((lag - 1) * k + 1:lag * k, :))
      end do
      do h = 1, horizon
         yhat = matmul(psi, d_pred(h, :))
         do lag = 1, min(h - 1, p)
            mu_fore = matmul(psi, d_pred(h - lag, :))
            yhat = yhat + matmul(pi_l(:, :, lag), y_pred(h - lag, :) - mu_fore)
         end do
         if (h <= p) then
            do lag = h, p
               mu_hist = matmul(psi, model%setup%x(n + h - lag, :))
               yhat = yhat + matmul(pi_l(:, :, lag), model%setup%y(n + h - lag, :) - mu_hist)
            end do
         end if
         call random_multivariate_normal(0.0_dp * yhat, sigma, innovation, info)
         if (info /= ssbvar_success) return
         y_pred(h, :) = yhat + innovation
      end do
      info = ssbvar_success
   end subroutine predictive_draw

   pure subroutine forecast_model(model, pi, use_median, growth_rate_idx, freq, result, info)
      type(bvar_model), intent(in) :: model !! Fitted BVAR object containing joint predictive draws.
      real(dp), intent(in), optional :: pi !! Central prediction-interval probability, strictly between zero and one; default 0.95.
      logical, intent(in), optional :: use_median !! Use posterior predictive medians instead of means for point forecasts.
      integer, intent(in), optional :: growth_rate_idx(:) !! One-based variables converted to trailing annual growth.
      integer, intent(in), optional :: freq !! Positive periods-per-year window used when growth-rate indices are supplied.
      type(forecast_result), intent(out) :: result !! Point, lower, and upper forecast matrices, each with shape (H,k).
      integer, intent(out) :: info !! Status code; zero indicates forecast summaries were computed.
      real(dp) :: interval
      logical :: med
      real(dp), allocatable :: draws(:, :, :)

      if (.not. model%has_fit) then
         info = ssbvar_not_ready
         return
      end if
      interval = 0.95_dp
      med = .false.
      if (present(pi)) interval = pi
      if (present(use_median)) med = use_median
      if (interval <= 0.0_dp .or. interval >= 1.0_dp) then
         info = ssbvar_invalid_input
         return
      end if
      allocate(draws(size(model%fit%y_pred, 1), size(model%fit%y_pred, 2), size(model%fit%y_pred, 3)))
      draws = model%fit%y_pred
      if (present(growth_rate_idx)) then
         if (.not. present(freq)) then
            info = ssbvar_invalid_input
            return
         end if
         call annualize_prediction_draws(model%data, draws, growth_rate_idx, freq, info)
         if (info /= ssbvar_success) return
      end if
      call summarize_prediction_draws(draws, interval, med, result)
      info = ssbvar_success
   end subroutine forecast_model

   subroutine conditional_forecast_model(model, condition_var, condition_horizon, condition_value, pi, &
      use_median, growth_rate_idx, freq, result, info)
      type(bvar_model), intent(in) :: model !! Homoscedastic fitted BVAR object whose posterior draws define conditional forecasts.
      integer, intent(in) :: condition_var(:) !! One-based variable index for each equality condition.
      integer, intent(in) :: condition_horizon(:) !! One-based forecast horizon for each equality condition, between 1 and H.
      real(dp), intent(in) :: condition_value(:) !! Target forecast value paired with each variable/horizon condition.
      real(dp), intent(in), optional :: pi !! Central conditional prediction-interval probability; default 0.95.
      logical, intent(in), optional :: use_median !! Use medians rather than means as conditional point forecasts.
      integer, intent(in), optional :: growth_rate_idx(:) !! One-based variables annualized by trailing sums.
      integer, intent(in), optional :: freq !! Positive trailing window size for annualized growth-rate variables.
      type(forecast_result), intent(out) :: result !! Conditional point, lower, and upper forecast matrices, each shape (H,k).
      integer, intent(out) :: info !! Status code; zero indicates conditional Gaussian shock draws succeeded.
      integer :: nd, hmax, k, p, v, s, n, row, h, j, lag, start_col
      real(dp) :: interval
      logical :: med
      real(dp), allocatable :: cond_draws(:, :, :), mean_path(:, :), pi_l(:, :, :), ma(:, :, :)
      real(dp), allocatable :: chol(:, :), rmat(:, :), rvec(:), rr(:, :), coeff(:), eta0(:), eta(:), z(:)
      real(dp), allocatable :: rz(:), correction(:), psi_weights(:, :, :), shock(:)

      if (.not. model%has_fit) then
         info = ssbvar_not_ready
         return
      end if
      if (.not. model%fit%homoscedastic .or. model%priors%sv) then
         info = ssbvar_unsupported
         return
      end if
      v = size(condition_var)
      if (size(condition_horizon) /= v .or. size(condition_value) /= v .or. v < 1) then
         info = ssbvar_invalid_input
         return
      end if
      nd = model%fit%n_draws
      hmax = model%fit%h
      k = model%setup%k
      p = model%setup%p
      s = k * hmax
      if (v > s) then
         info = ssbvar_invalid_input
         return
      end if
      do row = 1, v
         if (condition_var(row) < 1 .or. condition_var(row) > k) then
            info = ssbvar_invalid_input
            return
         end if
         if (condition_horizon(row) < 1 .or. condition_horizon(row) > hmax) then
            info = ssbvar_invalid_input
            return
         end if
      end do
      interval = 0.95_dp
      med = .false.
      if (present(pi)) interval = pi
      if (present(use_median)) med = use_median
      if (interval <= 0.0_dp .or. interval >= 1.0_dp) then
         info = ssbvar_invalid_input
         return
      end if

      allocate(cond_draws(nd, hmax, k), mean_path(hmax, k), pi_l(k, k, p))
      allocate(rmat(v, s), rvec(v), rr(v, v), coeff(v), eta0(s), eta(s), z(s))
      allocate(rz(v), correction(s), shock(k))
      do n = 1, nd
         call conditional_mean_path(model, model%fit%beta(:, :, n), model%fit%psi(:, :, n), mean_path)
         do lag = 1, p
            pi_l(:, :, lag) = transpose(model%fit%beta((lag - 1) * k + 1:lag * k, :, n))
         end do
         call var_psi_weights(pi_l, hmax - 1, psi_weights)
         call cholesky_factor(model%fit%sigma_u(:, :, n), chol, info)
         if (info /= 0) then
            info = ssbvar_linalg_failure
            return
         end if
         allocate(ma(k, k, 0:hmax - 1))
         do h = 0, hmax - 1
            ma(:, :, h) = matmul(psi_weights(:, :, h), chol)
         end do
         rmat = 0.0_dp
         do row = 1, v
            h = condition_horizon(row)
            do j = 1, h
               start_col = (j - 1) * k + 1
               rmat(row, start_col:start_col + k - 1) = ma(condition_var(row), :, h - j)
            end do
            rvec(row) = condition_value(row) - mean_path(h, condition_var(row))
         end do
         rr = matmul(rmat, transpose(rmat))
         call solve_spd(rr, rvec, coeff, info)
         if (info /= 0) then
            info = ssbvar_linalg_failure
            return
         end if
         eta0 = matmul(transpose(rmat), coeff)
         call random_standard_normal_vector(z)
         rz = matmul(rmat, z)
         call solve_spd(rr, rz, coeff, info)
         if (info /= 0) then
            info = ssbvar_linalg_failure
            return
         end if
         correction = matmul(transpose(rmat), coeff)
         eta = eta0 + z - correction
         do h = 1, hmax
            shock = 0.0_dp
            do j = 1, h
               shock = shock + matmul(ma(:, :, h - j), eta((j - 1) * k + 1:j * k))
            end do
            cond_draws(n, h, :) = mean_path(h, :) + shock
         end do
         deallocate(ma, chol, psi_weights)
      end do
      if (present(growth_rate_idx)) then
         if (.not. present(freq)) then
            info = ssbvar_invalid_input
            return
         end if
         call annualize_prediction_draws(model%data, cond_draws, growth_rate_idx, freq, info)
         if (info /= ssbvar_success) return
      end if
      call summarize_prediction_draws(cond_draws, interval, med, result)
      info = ssbvar_success
   end subroutine conditional_forecast_model

   subroutine irf_model(model, h, orthogonal, use_median, ci, growth_rate_idx, freq, result, info, t)
      type(bvar_model), intent(in) :: model !! Fitted BVAR object supplying posterior beta and covariance draws.
      integer, intent(in), optional :: h !! Maximum impulse-response horizon, inclusive of horizon zero; default 20.
      logical, intent(in), optional :: orthogonal !! True for Cholesky OIRFs; false for generalized IRFs.
      logical, intent(in), optional :: use_median !! Summarize posterior responses by medians rather than means.
      real(dp), intent(in), optional :: ci !! Central posterior credible-interval probability; default 0.95.
      integer, intent(in), optional :: growth_rate_idx(:) !! One-based response-variable indices to convert to trailing annual sums.
      integer, intent(in), optional :: freq !! Positive periods-per-year window required for annualized response variables.
      type(irf_result), intent(out) :: result !! Center, lower, and upper arrays with shape (k,k,H+1).
      integer, intent(out) :: info !! Status code; zero indicates posterior impulse responses were summarized.
      integer, intent(in), optional :: t !! Original-sample time selecting SV Sigma_u,t; default is the last observation.
      integer :: horizon, nd, k, p, n, lag, impulse, response, hh, t_est
      real(dp) :: interval, scale
      logical :: orth, med
      real(dp), allocatable :: draws(:, :, :, :), pi_l(:, :, :), weights(:, :, :), chol(:, :), impact(:, :)
      real(dp), allocatable :: path(:), annual(:), sigma(:, :)

      if (.not. model%has_fit) then
         info = ssbvar_not_ready
         return
      end if
      horizon = 20
      interval = 0.95_dp
      orth = .true.
      med = .true.
      if (present(h)) horizon = h
      if (present(ci)) interval = ci
      if (present(orthogonal)) orth = orthogonal
      if (present(use_median)) med = use_median
      if (horizon < 0 .or. interval <= 0.0_dp .or. interval >= 1.0_dp) then
         info = ssbvar_invalid_input
         return
      end if
      if (present(growth_rate_idx)) then
         if (.not. present(freq)) then
            info = ssbvar_invalid_input
            return
         end if
         if (freq <= 0) then
            info = ssbvar_invalid_input
            return
         end if
      end if

      nd = model%fit%n_draws
      k = model%setup%k
      p = model%setup%p
      t_est = model%setup%n
      if (.not. model%fit%homoscedastic) then
         if (present(t)) then
            if (t < p + 1 .or. t > p + model%setup%n) then
               info = ssbvar_invalid_input
               return
            end if
            t_est = t - p
         end if
      end if
      allocate(draws(nd, k, k, 0:horizon), pi_l(k, k, p), sigma(k, k))
      do n = 1, nd
         do lag = 1, p
            pi_l(:, :, lag) = transpose(model%fit%beta((lag - 1) * k + 1:lag * k, :, n))
         end do
         call var_psi_weights(pi_l, horizon, weights)
         if (model%fit%homoscedastic) then
            sigma = model%fit%sigma_u(:, :, n)
         else
            sigma = model%fit%sigma_u_time(:, :, t_est, n)
         end if
         if (orth) then
            call cholesky_factor(sigma, chol, info)
            if (info /= 0) then
               info = ssbvar_linalg_failure
               return
            end if
            do hh = 0, horizon
               impact = matmul(weights(:, :, hh), chol)
               draws(n, :, :, hh) = impact
            end do
            deallocate(chol)
         else
            do impulse = 1, k
               scale = sqrt(sigma(impulse, impulse))
               if (scale <= 0.0_dp) then
                  info = ssbvar_invalid_input
                  return
               end if
               do hh = 0, horizon
                  draws(n, :, impulse, hh) = matmul(weights(:, :, hh), sigma(:, impulse)) / scale
               end do
            end do
         end if
         deallocate(weights)
      end do

      if (present(growth_rate_idx)) then
         allocate(path(0:horizon), annual(0:horizon))
         do response = 1, k
            if (.not. integer_in_list(response, growth_rate_idx)) cycle
            do impulse = 1, k
               do n = 1, nd
                  path = draws(n, response, impulse, :)
                  do hh = 0, horizon
                     annual(hh) = sum(path(max(0, hh - freq + 1):hh))
                  end do
                  draws(n, response, impulse, :) = annual
               end do
            end do
         end do
      end if
      call summarize_irf_draws(draws, interval, med, result)
      info = ssbvar_success
   end subroutine irf_model

   pure subroutine bvar_summary(model, use_median, result, info, t)
      type(bvar_model), intent(in) :: model !! Fitted BVAR object whose posterior parameter summaries are requested.
      logical, intent(in), optional :: use_median !! Return elementwise posterior medians when true, otherwise means.
      type(bvar_summary_result), intent(out) :: result !! Posterior summaries for the fitted volatility specification.
      integer, intent(out) :: info !! Status code; zero indicates summary matrices were constructed successfully.
      integer, intent(in), optional :: t !! Original-sample time selecting SV Sigma_u,t; default is the last observation.
      integer :: t_est
      logical :: med

      if (.not. model%has_fit) then
         info = ssbvar_not_ready
         return
      end if
      med = .false.
      if (present(use_median)) med = use_median
      if (med) then
         result%beta = model%fit%beta_median
         result%psi = model%fit%psi_median
      else
         result%beta = model%fit%beta_mean
         result%psi = model%fit%psi_mean
      end if

      if (model%fit%homoscedastic) then
         if (med) then
            result%sigma_u = model%fit%sigma_median
         else
            result%sigma_u = model%fit%sigma_mean
         end if
         info = ssbvar_success
         return
      end if

      t_est = model%setup%n
      if (present(t)) then
         if (t < model%setup%p + 1 .or. t > model%setup%p + model%setup%n) then
            info = ssbvar_invalid_input
            return
         end if
         t_est = t - model%setup%p
      end if
      if (med) then
         result%sigma_u = model%fit%sigma_time_median(:, :, t_est)
      else
         result%sigma_u = model%fit%sigma_time_mean(:, :, t_est)
      end if
      if (med) then
         result%a = model%fit%a_median
      else
         result%a = model%fit%a_mean
      end if
      if (trim(model%fit%sv_type) == "RW") then
         if (med) then
            result%phi = model%fit%phi_median
         else
            result%phi = model%fit%phi_mean
         end if
      else if (trim(model%fit%sv_type) == "AR1") then
         if (med) then
            result%gamma_0 = model%fit%gamma_0_median
            result%gamma_1 = model%fit%gamma_1_median
            result%phi_cov = model%fit%phi_cov_median
         else
            result%gamma_0 = model%fit%gamma_0_mean
            result%gamma_1 = model%fit%gamma_1_mean
            result%phi_cov = model%fit%phi_cov_mean
         end if
      end if
      info = ssbvar_success
   end subroutine bvar_summary

   pure subroutine build_future_d(model, horizon, d_pred, future_d, info)
      type(bvar_model), intent(in) :: model !! Set-up BVAR object providing the observed deterministic matrix and specification.
      integer, intent(in) :: horizon !! Positive number of future deterministic rows to construct or validate.
      real(dp), intent(in), optional :: d_pred(:, :) !! Optional caller-supplied future deterministic values with shape (H,q).
      real(dp), allocatable, intent(out) :: future_d(:, :) !! Allocated future deterministic matrix with shape (H,q).
      integer, intent(out) :: info !! Status code; zero indicates successful construction or validation.
      integer :: q, h
      real(dp), allocatable :: last_row(:)

      q = model%setup%q
      allocate(future_d(horizon, q), last_row(q))
      if (present(d_pred)) then
         if (size(d_pred, 1) /= horizon .or. size(d_pred, 2) /= q) then
            info = ssbvar_invalid_input
            return
         end if
         future_d = d_pred
      else
         last_row = model%setup%dt(size(model%setup%dt, 1), :)
         do h = 1, horizon
            future_d(h, :) = last_row
         end do
         if (trim(model%setup%deterministic) == "constant_and_trend") then
            do h = 1, horizon
               future_d(h, 2) = last_row(2) + real(h, dp)
            end do
         end if
      end if
      info = ssbvar_success
   end subroutine build_future_d

   pure subroutine conditional_mean_path(model, beta, psi, path)
      type(bvar_model), intent(in) :: model !! Fitted model and observed sample supplying pre-forecast lags.
      real(dp), intent(in) :: beta(:, :) !! Autoregressive coefficient draw with shape (k*p,k).
      real(dp), intent(in) :: psi(:, :) !! Steady-state parameter draw with shape (k,q).
      real(dp), intent(out) :: path(:, :) !! Shock-free forecast path with shape (H,k), where H is the fitted forecast horizon.
      integer :: horizon, k, p, n, h, lag
      real(dp), allocatable :: pi_l(:, :, :), mu(:)

      horizon = model%fit%h
      k = model%setup%k
      p = model%setup%p
      n = model%setup%n
      allocate(pi_l(k, k, p), mu(k))
      do lag = 1, p
         pi_l(:, :, lag) = transpose(beta((lag - 1) * k + 1:lag * k, :))
      end do
      do h = 1, horizon
         path(h, :) = matmul(psi, model%fit%d_pred(h, :))
         do lag = 1, min(h - 1, p)
            mu = matmul(psi, model%fit%d_pred(h - lag, :))
            path(h, :) = path(h, :) + matmul(pi_l(:, :, lag), path(h - lag, :) - mu)
         end do
         if (h <= p) then
            do lag = h, p
               mu = matmul(psi, model%setup%x(n + h - lag, :))
               path(h, :) = path(h, :) + &
                  matmul(pi_l(:, :, lag), model%setup%y(n + h - lag, :) - mu)
            end do
         end if
      end do
   end subroutine conditional_mean_path

   pure subroutine compute_fit_summaries(model)
      type(bvar_model), intent(inout) :: model !! BVAR object whose retained draws are reduced to elementwise means and medians.
      integer :: kp, k, q, i, j, t
      real(dp), allocatable :: work(:)

      kp = size(model%fit%beta, 1)
      k = size(model%fit%beta, 2)
      q = size(model%fit%psi, 2)
      allocate(model%fit%beta_mean(kp, k), model%fit%beta_median(kp, k))
      allocate(model%fit%psi_mean(k, q), model%fit%psi_median(k, q))
      allocate(work(model%fit%n_draws))
      do j = 1, k
         do i = 1, kp
            work = model%fit%beta(i, j, :)
            model%fit%beta_mean(i, j) = sum(work) / real(size(work), dp)
            model%fit%beta_median(i, j) = r_quantile_type7(work, 0.5_dp)
         end do
      end do
      do j = 1, q
         do i = 1, k
            work = model%fit%psi(i, j, :)
            model%fit%psi_mean(i, j) = sum(work) / real(size(work), dp)
            model%fit%psi_median(i, j) = r_quantile_type7(work, 0.5_dp)
         end do
      end do

      if (model%fit%homoscedastic) then
         allocate(model%fit%sigma_mean(k, k), model%fit%sigma_median(k, k))
         do j = 1, k
            do i = 1, k
               work = model%fit%sigma_u(i, j, :)
               model%fit%sigma_mean(i, j) = sum(work) / real(size(work), dp)
               model%fit%sigma_median(i, j) = r_quantile_type7(work, 0.5_dp)
            end do
         end do
         return
      end if

      allocate(model%fit%sigma_time_mean(k, k, model%setup%n))
      allocate(model%fit%sigma_time_median(k, k, model%setup%n))
      do t = 1, model%setup%n
         do j = 1, k
            do i = 1, k
               work = model%fit%sigma_u_time(i, j, t, :)
               model%fit%sigma_time_mean(i, j, t) = sum(work) / real(size(work), dp)
               model%fit%sigma_time_median(i, j, t) = r_quantile_type7(work, 0.5_dp)
            end do
         end do
      end do

      allocate(model%fit%a_mean(k, k), model%fit%a_median(k, k))
      do j = 1, k
         do i = 1, k
            work = model%fit%a(i, j, :)
            model%fit%a_mean(i, j) = sum(work) / real(size(work), dp)
            model%fit%a_median(i, j) = r_quantile_type7(work, 0.5_dp)
         end do
      end do
      if (trim(model%fit%sv_type) == "RW") then
         allocate(model%fit%phi_mean(k), model%fit%phi_median(k))
         do i = 1, k
            work = model%fit%phi(i, :)
            model%fit%phi_mean(i) = sum(work) / real(size(work), dp)
            model%fit%phi_median(i) = r_quantile_type7(work, 0.5_dp)
         end do
      else if (trim(model%fit%sv_type) == "AR1") then
         allocate(model%fit%gamma_0_mean(k), model%fit%gamma_0_median(k))
         allocate(model%fit%gamma_1_mean(k), model%fit%gamma_1_median(k))
         allocate(model%fit%phi_cov_mean(k, k), model%fit%phi_cov_median(k, k))
         do i = 1, k
            work = model%fit%gamma_0(i, :)
            model%fit%gamma_0_mean(i) = sum(work) / real(size(work), dp)
            model%fit%gamma_0_median(i) = r_quantile_type7(work, 0.5_dp)
            work = model%fit%gamma_1(i, :)
            model%fit%gamma_1_mean(i) = sum(work) / real(size(work), dp)
            model%fit%gamma_1_median(i) = r_quantile_type7(work, 0.5_dp)
         end do
         do j = 1, k
            do i = 1, k
               work = model%fit%phi_cov(i, j, :)
               model%fit%phi_cov_mean(i, j) = sum(work) / real(size(work), dp)
               model%fit%phi_cov_median(i, j) = r_quantile_type7(work, 0.5_dp)
            end do
         end do
      end if
   end subroutine compute_fit_summaries

   pure subroutine summarize_prediction_draws(draws, interval, use_median, result)
      real(dp), intent(in) :: draws(:, :, :) !! Predictive draws with shape (draw,H,k).
      real(dp), intent(in) :: interval !! Central interval probability in the open unit interval.
      logical, intent(in) :: use_median !! Select median rather than mean for the point summary.
      type(forecast_result), intent(out) :: result !! Allocated point, lower, and upper forecast matrices.
      integer :: h, k
      real(dp) :: alpha
      real(dp), allocatable :: work(:)

      allocate(result%forecast(size(draws, 2), size(draws, 3)))
      allocate(result%lower(size(draws, 2), size(draws, 3)))
      allocate(result%upper(size(draws, 2), size(draws, 3)))
      allocate(work(size(draws, 1)))
      alpha = 1.0_dp - interval
      do k = 1, size(draws, 3)
         do h = 1, size(draws, 2)
            work = draws(:, h, k)
            if (use_median) then
               result%forecast(h, k) = r_quantile_type7(work, 0.5_dp)
            else
               result%forecast(h, k) = sum(work) / real(size(work), dp)
            end if
            result%lower(h, k) = r_quantile_type7(work, alpha / 2.0_dp)
            result%upper(h, k) = r_quantile_type7(work, 1.0_dp - alpha / 2.0_dp)
         end do
      end do
   end subroutine summarize_prediction_draws

   pure subroutine summarize_irf_draws(draws, interval, use_median, result)
      real(dp), intent(in) :: draws(:, :, :, :) !! Posterior IRF draws with shape (draw,response,impulse,horizon).
      real(dp), intent(in) :: interval !! Central posterior credible-interval probability.
      logical, intent(in) :: use_median !! Select median rather than mean for the center response.
      type(irf_result), intent(out) :: result !! Allocated center, lower, and upper arrays with the non-draw dimensions.
      integer :: response, impulse, h, h_lower, h_upper
      real(dp) :: alpha
      real(dp), allocatable :: work(:)

      h_lower = lbound(draws, 4)
      h_upper = ubound(draws, 4)
      allocate(result%center(size(draws, 2), size(draws, 3), h_lower:h_upper))
      allocate(result%lower(size(draws, 2), size(draws, 3), h_lower:h_upper))
      allocate(result%upper(size(draws, 2), size(draws, 3), h_lower:h_upper))
      allocate(work(size(draws, 1)))
      alpha = 1.0_dp - interval
      do impulse = 1, size(draws, 3)
         do response = 1, size(draws, 2)
            do h = h_lower, h_upper
               work = draws(:, response, impulse, h)
               if (use_median) then
                  result%center(response, impulse, h) = r_quantile_type7(work, 0.5_dp)
               else
                  result%center(response, impulse, h) = sum(work) / real(size(work), dp)
               end if
               result%lower(response, impulse, h) = r_quantile_type7(work, alpha / 2.0_dp)
               result%upper(response, impulse, h) = r_quantile_type7(work, 1.0_dp - alpha / 2.0_dp)
            end do
         end do
      end do
   end subroutine summarize_irf_draws

   pure subroutine annualize_prediction_draws(history, draws, indices, freq, info)
      real(dp), intent(in) :: history(:, :) !! Observed time-by-variable data supplying pre-forecast growth rates.
      real(dp), intent(inout) :: draws(:, :, :) !! Predictive draws (draw,H,k) modified in place for selected variables.
      integer, intent(in) :: indices(:) !! One-based variables to replace by trailing sums over `freq` periods.
      integer, intent(in) :: freq !! Positive trailing-window length, normally 4 for quarterly or 12 for monthly data.
      integer, intent(out) :: info !! Status code; zero indicates all requested annualizations were performed.
      integer :: n_draws, horizon, k, var, draw, h, r, n_hist, offset
      real(dp), allocatable :: extended(:)

      if (freq <= 0) then
         info = ssbvar_invalid_input
         return
      end if
      n_hist = size(history, 1)
      if (n_hist < freq - 1) then
         info = ssbvar_invalid_input
         return
      end if
      n_draws = size(draws, 1)
      horizon = size(draws, 2)
      k = size(draws, 3)
      allocate(extended(freq - 1 + horizon))
      do r = 1, size(indices)
         var = indices(r)
         if (var < 1 .or. var > k) then
            info = ssbvar_invalid_input
            return
         end if
         do draw = 1, n_draws
            if (freq > 1) then
               extended(1:freq - 1) = history(n_hist - freq + 2:n_hist, var)
            end if
            extended(freq:freq - 1 + horizon) = draws(draw, :, var)
            do h = 1, horizon
               offset = h
               draws(draw, h, var) = sum(extended(offset:offset + freq - 1))
            end do
         end do
      end do
      info = ssbvar_success
   end subroutine annualize_prediction_draws

   pure logical function integer_in_list(value, values) result(found)
      integer, intent(in) :: value !! Integer whose membership is tested.
      integer, intent(in) :: values(:) !! Candidate integer list; duplicates are harmless.
      found = any(values == value)
   end function integer_in_list

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n !! Nonnegative matrix order.
      real(dp) :: a(n, n)
      integer :: i

      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

   subroutine seed_rng(seed)
      integer, intent(in) :: seed !! Deterministic scalar seed expanded to the compiler-required intrinsic RNG state length.
      integer :: n, i
      integer, allocatable :: put(:)

      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729 * i + 8191 * i * i, huge(1) - 1)
         if (put(i) == 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_rng

   subroutine random_standard_normal_vector(z)
      real(dp), intent(out) :: z(:) !! Independent standard-normal variates filled using Box-Muller transforms.
      integer :: i
      real(dp) :: u1, u2, radius, angle
      real(dp), parameter :: two_pi = 6.283185307179586_dp

      i = 1
      do while (i <= size(z))
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         radius = sqrt(-2.0_dp * log(u1))
         angle = two_pi * u2
         z(i) = radius * cos(angle)
         if (i + 1 <= size(z)) z(i + 1) = radius * sin(angle)
         i = i + 2
      end do
   end subroutine random_standard_normal_vector

   subroutine random_multivariate_normal(mean_value, covariance, draw, info)
      real(dp), intent(in) :: mean_value(:) !! Mean vector of the multivariate normal distribution.
      real(dp), intent(in) :: covariance(:, :) !! Symmetric positive-definite covariance matrix conforming with the mean.
      real(dp), intent(out) :: draw(:) !! Generated multivariate normal vector with the same length as the mean.
      integer, intent(out) :: info !! Status code; zero indicates Cholesky factorization and sampling succeeded.
      real(dp), allocatable :: factor(:, :), z(:)

      if (size(covariance, 1) /= size(mean_value) .or. size(covariance, 2) /= size(mean_value) .or. &
          size(draw) /= size(mean_value)) then
         info = ssbvar_invalid_input
         return
      end if
      call cholesky_factor(covariance, factor, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      allocate(z(size(mean_value)))
      call random_standard_normal_vector(z)
      draw = mean_value + matmul(factor, z)
      info = ssbvar_success
   end subroutine random_multivariate_normal

   recursive subroutine random_gamma(shape, value)
      real(dp), intent(in) :: shape !! Positive gamma shape parameter for a unit-scale variate.
      real(dp), intent(out) :: value !! Generated unit-scale gamma variate.
      real(dp) :: d, c, x, v, u, small
      real(dp) :: z(1)

      if (shape < 1.0_dp) then
         call random_gamma(shape + 1.0_dp, small)
         call random_number(u)
         u = max(u, tiny(1.0_dp))
         value = small * u ** (1.0_dp / shape)
         return
      end if
      d = shape - 1.0_dp / 3.0_dp
      c = 1.0_dp / sqrt(9.0_dp * d)
      do
         call random_standard_normal_vector(z)
         x = z(1)
         v = 1.0_dp + c * x
         if (v <= 0.0_dp) cycle
         v = v * v * v
         call random_number(u)
         if (u < 1.0_dp - 0.0331_dp * x ** 4) exit
         if (log(max(u, tiny(1.0_dp))) < 0.5_dp * x * x + d * (1.0_dp - v + log(v))) exit
      end do
      value = d * v
   end subroutine random_gamma

   subroutine random_chisquare(df, value)
      real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.
      real(dp), intent(out) :: value !! Generated chi-square variate.
      real(dp) :: gamma_value

      call random_gamma(0.5_dp * df, gamma_value)
      value = 2.0_dp * gamma_value
   end subroutine random_chisquare

   subroutine random_inverse_wishart(df, scale, draw, info)
      integer, intent(in) :: df !! Inverse-Wishart degrees of freedom, greater than matrix order minus one.
      real(dp), intent(in) :: scale(:, :) !! Symmetric positive-definite inverse-Wishart scale matrix.
      real(dp), intent(out) :: draw(:, :) !! Generated inverse-Wishart covariance matrix conforming with `scale`.
      integer, intent(out) :: info !! Status code; zero indicates the Bartlett construction and inversion succeeded.
      integer :: k, i, j
      real(dp), allocatable :: scale_inv(:, :), factor(:, :), bartlett(:, :), work(:, :), wishart(:, :)
      real(dp) :: chi
      real(dp) :: z(1)

      k = size(scale, 1)
      if (size(scale, 2) /= k .or. size(draw, 1) /= k .or. size(draw, 2) /= k .or. df <= k - 1) then
         info = ssbvar_invalid_input
         return
      end if
      call inverse_matrix(scale, scale_inv, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      call cholesky_factor(scale_inv, factor, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      allocate(bartlett(k, k), work(k, k), wishart(k, k))
      bartlett = 0.0_dp
      do i = 1, k
         call random_chisquare(real(df - i + 1, dp), chi)
         bartlett(i, i) = sqrt(chi)
         do j = 1, i - 1
            call random_standard_normal_vector(z)
            bartlett(i, j) = z(1)
         end do
      end do
      work = matmul(factor, bartlett)
      wishart = matmul(work, transpose(work))
      call inverse_matrix(wishart, scale_inv, info)
      if (info /= 0) then
         info = ssbvar_linalg_failure
         return
      end if
      draw = 0.5_dp * (scale_inv + transpose(scale_inv))
      info = ssbvar_success
   end subroutine random_inverse_wishart

end module steadystatebvar_core
