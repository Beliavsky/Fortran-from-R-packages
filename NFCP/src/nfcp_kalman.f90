module nfcp_kalman
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use nfcp_types, only : dp, nfcp_model_t, nfcp_filter_result_t, nfcp_ok, &
                         nfcp_invalid_input, nfcp_singular, nfcp_nonfinite
  use nfcp_math, only : nfcp_pi, nfcp_a_t, nfcp_covariance, nfcp_seasonality, &
                        inverse_spd, logdet_spd
  implicit none
  private
  public :: nfcp_kalman_filter

contains

  subroutine nfcp_kalman_filter(model, log_futures, dt, futures_ttm, result, &
                                initial_state, me_ttm, seasonal_trend, n_parameters)
    type(nfcp_model_t), intent(in) :: model
    real(dp), intent(in) :: log_futures(:,:)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: futures_ttm(:,:)
    type(nfcp_filter_result_t), intent(out) :: result
    real(dp), intent(in), optional :: initial_state(:)
    real(dp), intent(in), optional :: me_ttm(:)
    real(dp), intent(in), optional :: seasonal_trend
    integer, intent(in), optional :: n_parameters

    integer :: nt, nc, nf, t, i, j, m, status, npar, n_scalar_obs
    integer, allocatable :: obs(:)
    real(dp) :: trend, logdet, quad, ll_inc, nanv
    real(dp), allocatable :: x(:), xp(:), p(:,:), pp(:,:), q(:,:), g(:,:), c(:)
    real(dp), allocatable :: z(:,:), h(:,:), f(:,:), finv(:,:), k_gain(:,:)
    real(dp), allocatable :: y(:), d(:), innovation(:), fitted(:), temp(:,:), ident(:,:)
    character(len=:), allocatable :: validation_message

    nanv = ieee_value(0.0_dp, ieee_quiet_nan)
    result%status = nfcp_invalid_input
    result%message = 'invalid input'

    if (.not. model%valid(validation_message)) then
      result%message = validation_message
      return
    end if
    nt = size(log_futures,1)
    nc = size(log_futures,2)
    nf = model%n_factors
    if (nt < 1 .or. nc < 1 .or. dt <= 0.0_dp) then
      result%message = 'log_futures must be nonempty and dt positive'
      return
    end if
    if (any(shape(futures_ttm) /= shape(log_futures))) then
      result%message = 'futures_ttm must have the same shape as log_futures'
      return
    end if
    if (present(initial_state)) then
      if (size(initial_state) /= nf .or. .not. all(ieee_is_finite(initial_state))) then
        result%message = 'initial_state has invalid size or values'
        return
      end if
    end if
    if (model%n_me > 1 .and. model%n_me < nc) then
      if (.not. present(me_ttm)) then
        result%message = 'me_ttm is required for maturity-grouped measurement errors'
        return
      end if
      if (size(me_ttm) /= model%n_me .or. any(me_ttm(2:) <= me_ttm(:model%n_me-1))) then
        result%message = 'me_ttm must be strictly increasing with length n_me'
        return
      end if
    end if

    trend = 0.0_dp
    if (present(seasonal_trend)) trend = seasonal_trend
    if (trend < 0.0_dp .or. trend > 1.0_dp) then
      result%message = 'seasonal_trend must be between zero and one'
      return
    end if

    allocate(x(nf), xp(nf), p(nf,nf), pp(nf,nf), q(nf,nf), g(nf,nf), c(nf))
    allocate(ident(nf,nf), obs(nc))
    allocate(result%final_state(nf), result%final_covariance(nf,nf))
    allocate(result%states(nt,nf), result%state_variances(nt,nf))
    allocate(result%fitted_log_futures(nt,nc), result%residuals(nt,nc))
    allocate(result%log_likelihood_path(nt))
    result%fitted_log_futures = nanv
    result%residuals = nanv
    result%states = nanv
    result%state_variances = nanv
    result%log_likelihood_path = 0.0_dp

    x = 0.0_dp
    if (present(initial_state)) then
      x = initial_state
    else if (model%gbm) then
      do i = 1, nc
        if (ieee_is_finite(log_futures(1,i))) then
          x(1) = log_futures(1,i)
          exit
        end if
      end do
    end if
    p = 0.0_dp
    ident = 0.0_dp
    do i = 1, nf
      p(i,i) = 100.0_dp
      ident(i,i) = 1.0_dp
    end do
    call nfcp_covariance(model, dt, q)
    g = 0.0_dp
    do i = 1, nf
      g(i,i) = exp(-model%kappa(i)*dt)
    end do
    c = 0.0_dp
    if (model%gbm) c(1) = model%mu*dt

    result%log_likelihood = 0.0_dp
    n_scalar_obs = 0
    do t = 1, nt
      xp = c + matmul(g,x)
      pp = matmul(g,matmul(p,transpose(g))) + q
      pp = 0.5_dp*(pp + transpose(pp))

      m = 0
      do i = 1, nc
        if (ieee_is_finite(log_futures(t,i)) .and. ieee_is_finite(futures_ttm(t,i))) then
          m = m + 1
          obs(m) = i
        end if
      end do

      if (m > 0) then
        allocate(z(m,nf), h(m,m), f(m,m), finv(m,m), k_gain(nf,m))
        allocate(y(m), d(m), innovation(m), fitted(m), temp(nf,nf))
        z = 0.0_dp
        h = 0.0_dp
        do i = 1, m
          j = obs(i)
          y(i) = log_futures(t,j)
          d(i) = model%equilibrium + nfcp_seasonality(model, futures_ttm(t,j) + trend + real(t-1,dp)*dt) + &
                 nfcp_a_t(model, futures_ttm(t,j))
          z(i,:) = exp(-model%kappa*futures_ttm(t,j))
          h(i,i) = measurement_variance(model, j, futures_ttm(t,j), nc, me_ttm)
        end do
        f = matmul(z,matmul(pp,transpose(z))) + h
        f = 0.5_dp*(f + transpose(f))
        call inverse_spd(f, finv, status)
        if (status /= nfcp_ok) then
          result%status = nfcp_singular
          result%message = 'singular innovation covariance in Kalman filter'
          return
        end if
        call logdet_spd(f, logdet, status)
        if (status /= nfcp_ok .or. .not. ieee_is_finite(logdet)) then
          result%status = nfcp_singular
          result%message = 'non-positive innovation covariance in Kalman filter'
          return
        end if
        fitted = d + matmul(z,xp)
        innovation = y-fitted
        quad = dot_product(innovation,matmul(finv,innovation))
        if (.not. ieee_is_finite(quad)) then
          result%status = nfcp_nonfinite
          result%message = 'non-finite Kalman likelihood contribution'
          return
        end if
        ll_inc = -0.5_dp*(real(m,dp)*log(2.0_dp*nfcp_pi) + logdet + quad)
        result%log_likelihood = result%log_likelihood + ll_inc
        n_scalar_obs = n_scalar_obs + m

        k_gain = matmul(pp,matmul(transpose(z),finv))
        x = xp + matmul(k_gain,innovation)
        temp = ident - matmul(k_gain,z)
        p = matmul(temp,matmul(pp,transpose(temp))) + matmul(k_gain,matmul(h,transpose(k_gain)))
        p = 0.5_dp*(p + transpose(p))

        fitted = d + matmul(z,x)
        do i = 1, m
          j = obs(i)
          result%fitted_log_futures(t,j) = fitted(i)
          result%residuals(t,j) = fitted(i)-y(i)
        end do
        deallocate(z,h,f,finv,k_gain,y,d,innovation,fitted,temp)
      else
        x = xp
        p = pp
      end if

      result%states(t,:) = x
      do i = 1, nf
        result%state_variances(t,i) = p(i,i)
      end do
      result%log_likelihood_path(t) = result%log_likelihood
    end do

    result%final_state = x
    result%final_covariance = p
    npar = 0
    if (present(n_parameters)) npar = n_parameters
    result%aic = 2.0_dp*real(npar,dp)-2.0_dp*result%log_likelihood
    if (n_scalar_obs > 0) then
      result%bic = real(npar,dp)*log(real(n_scalar_obs,dp))-2.0_dp*result%log_likelihood
    end if
    result%status = nfcp_ok
    result%message = 'ok'
  end subroutine nfcp_kalman_filter

  real(dp) function measurement_variance(model, contract, maturity, n_contracts, me_ttm) result(v)
    type(nfcp_model_t), intent(in) :: model
    integer, intent(in) :: contract, n_contracts
    real(dp), intent(in) :: maturity
    real(dp), intent(in), optional :: me_ttm(:)
    integer :: k

    if (model%n_me <= 0) then
      v = 0.0_dp
    else if (model%n_me == 1) then
      v = model%measurement_error(1)**2
    else if (model%n_me == n_contracts) then
      v = model%measurement_error(contract)**2
    else
      k = model%n_me
      if (present(me_ttm)) then
        do k = 1, model%n_me
          if (maturity < me_ttm(k)) exit
        end do
        k = min(k,model%n_me)
      end if
      v = model%measurement_error(k)**2
    end if
    if (v < 1.01e-10_dp) v = 0.0_dp
  end function measurement_variance

end module nfcp_kalman
