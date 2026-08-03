! SPDX-License-Identifier: GPL-3.0-only
module rumidas_garch_midas
  use rumidas_kinds, only: dp
  use rumidas_status
  use rumidas_types
  use rumidas_weights, only: midas_weighted_component
  use rumidas_statistics, only: sample_mean
  implicit none
  private

  real(dp), parameter :: LOG_TWO_PI = log(2.0_dp * acos(-1.0_dp))

  public :: garch_midas_parameter_count, garch_midas_evaluate
  public :: gm_loglik, gm_cond_vol, gm_long_run_vol
  public :: gm_loglik_no_skew, gm_cond_vol_no_skew, gm_long_run_vol_no_skew
  public :: gm_x_loglik, gm_x_cond_vol, gm_x_long_run_vol
  public :: gm_x_loglik_no_skew, gm_x_cond_vol_no_skew, gm_x_long_run_vol_no_skew
  public :: gm_2m_loglik, gm_2m_cond_vol, gm_2m_long_run_vol
  public :: gm_2m_loglik_no_skew, gm_2m_cond_vol_no_skew, gm_2m_long_run_vol_no_skew
  public :: dagm_loglik, dagm_cond_vol, dagm_long_run_vol
  public :: dagm_loglik_no_skew, dagm_cond_vol_no_skew, dagm_long_run_vol_no_skew
  public :: dagm_x_loglik, dagm_x_cond_vol, dagm_x_long_run_vol
  public :: dagm_x_loglik_no_skew, dagm_x_cond_vol_no_skew, dagm_x_long_run_vol_no_skew
  public :: dagm_2m_loglik, dagm_2m_cond_vol, dagm_2m_long_run
  public :: dagm_2m_loglik_no_skew, dagm_2m_cond_vol_no_skew, dagm_2m_long_run_no_skew

contains

  pure integer function garch_midas_parameter_count(spec) result(number)
    type(garch_midas_spec), intent(in) :: spec
    number = 2
    if (spec%skew) number = number + 1
    if (spec%model == RUMIDAS_GMX .or. spec%model == RUMIDAS_DAGMX) number = number + 1
    number = number + 1
    select case (spec%model)
    case (RUMIDAS_GM, RUMIDAS_GMX)
      number = number + 2
    case (RUMIDAS_GM2M)
      number = number + 4
    case (RUMIDAS_DAGM, RUMIDAS_DAGMX)
      number = number + 4
    case (RUMIDAS_DAGM2M)
      number = number + 8
    case default
      number = 0
      return
    end select
    if (spec%distribution == RUMIDAS_STUDENT_T) number = number + 1
  end function garch_midas_parameter_count

  subroutine garch_midas_evaluate(param, spec, daily_ret, mv_m_1, loglik, conditional_volatility, &
      long_run_volatility, short_run_variance, status, mv_m_2, x_variable)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :)
    type(garch_midas_spec), intent(in) :: spec
    real(dp), allocatable, intent(out) :: loglik(:), conditional_volatility(:)
    real(dp), allocatable, intent(out) :: long_run_volatility(:), short_run_variance(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: mv_m_2(:, :), x_variable(:)

    real(dp), allocatable :: tau(:), c1(:), c2(:), c3(:), c4(:), step(:)
    real(dp) :: alpha, beta, gamma_value, x_coefficient, intercept, m, degrees
    real(dp) :: theta(4), w2(4), variance, residual, mean_return
    integer :: n, i, index, local_status, expected
    logical :: has_x, two_midas, double_asymmetric

    status = RUMIDAS_SUCCESS
    n = size(daily_ret)
    expected = garch_midas_parameter_count(spec)
    if (n <= 0 .or. expected <= 0 .or. size(param) /= expected) then
      call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
      status = RUMIDAS_INVALID_INPUT
      return
    end if
    if (size(mv_m_1, 1) /= spec%k1 + 1 .or. size(mv_m_1, 2) /= n) then
      call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
      status = RUMIDAS_DIMENSION_ERROR
      return
    end if

    has_x = spec%model == RUMIDAS_GMX .or. spec%model == RUMIDAS_DAGMX
    two_midas = spec%model == RUMIDAS_GM2M .or. spec%model == RUMIDAS_DAGM2M
    double_asymmetric = spec%model == RUMIDAS_DAGM .or. spec%model == RUMIDAS_DAGMX .or. &
      spec%model == RUMIDAS_DAGM2M
    if (has_x) then
      if (.not. present(x_variable)) then
        call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
        status = RUMIDAS_INVALID_INPUT
        return
      end if
      if (size(x_variable) /= n) then
        call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
        status = RUMIDAS_DIMENSION_ERROR
        return
      end if
    end if
    if (two_midas) then
      if (.not. present(mv_m_2)) then
        call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
        status = RUMIDAS_INVALID_INPUT
        return
      end if
      if (size(mv_m_2, 1) /= spec%k2 + 1 .or. size(mv_m_2, 2) /= n) then
        call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
        status = RUMIDAS_DIMENSION_ERROR
        return
      end if
    end if

    alpha = param(1)
    beta = param(2)
    index = 3
    gamma_value = 0.0_dp
    if (spec%skew) then
      gamma_value = param(index)
      index = index + 1
    end if
    x_coefficient = 0.0_dp
    if (has_x) then
      x_coefficient = param(index)
      index = index + 1
    end if
    m = param(index)
    index = index + 1
    theta = 0.0_dp
    w2 = 0.0_dp
    select case (spec%model)
    case (RUMIDAS_GM, RUMIDAS_GMX)
      theta(1) = param(index); w2(1) = param(index + 1); index = index + 2
    case (RUMIDAS_GM2M)
      theta(1) = param(index); w2(1) = param(index + 1)
      theta(2) = param(index + 2); w2(2) = param(index + 3); index = index + 4
    case (RUMIDAS_DAGM, RUMIDAS_DAGMX)
      theta(1) = param(index); w2(1) = param(index + 1)
      theta(2) = param(index + 2); w2(2) = param(index + 3); index = index + 4
    case (RUMIDAS_DAGM2M)
      theta(1) = param(index); w2(1) = param(index + 1)
      theta(2) = param(index + 2); w2(2) = param(index + 3)
      theta(3) = param(index + 4); w2(3) = param(index + 5)
      theta(4) = param(index + 6); w2(4) = param(index + 7); index = index + 8
    end select
    degrees = 0.0_dp
    if (spec%distribution == RUMIDAS_STUDENT_T) degrees = param(index)

    if (.not. valid_garch_parameters(alpha, beta, gamma_value, spec%skew, degrees, spec%distribution)) then
      call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
      status = RUMIDAS_INVALID_PARAMETER
      return
    end if
    if (.not. valid_lag_parameters(w2, spec)) then
      call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
      status = RUMIDAS_INVALID_PARAMETER
      return
    end if

    allocate(tau(n), c1(n), c2(n), c3(n), c4(n), step(n))
    c1 = 0.0_dp; c2 = 0.0_dp; c3 = 0.0_dp; c4 = 0.0_dp
    if (.not. double_asymmetric) then
      call midas_weighted_component(mv_m_1, spec%k1, w2(1), spec%lag_function, c1, local_status)
      if (local_status /= RUMIDAS_SUCCESS) then
        status = local_status
        call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
        return
      end if
      if (two_midas) then
        call midas_weighted_component(mv_m_2, spec%k2, w2(2), spec%lag_function, c2, local_status)
        if (local_status /= RUMIDAS_SUCCESS) then
          status = local_status
          call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
          return
        end if
      end if
    else
      call midas_weighted_component(mv_m_1, spec%k1, w2(1), spec%lag_function, c1, local_status, 1)
      call midas_weighted_component(mv_m_1, spec%k1, w2(2), spec%lag_function, c2, local_status, -1)
      if (local_status /= RUMIDAS_SUCCESS) then
        status = local_status
        call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
        return
      end if
      if (two_midas) then
        call midas_weighted_component(mv_m_2, spec%k2, w2(3), spec%lag_function, c3, local_status, 1)
        call midas_weighted_component(mv_m_2, spec%k2, w2(4), spec%lag_function, c4, local_status, -1)
        if (local_status /= RUMIDAS_SUCCESS) then
          status = local_status
          call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
          return
        end if
      end if
    end if

    select case (spec%model)
    case (RUMIDAS_GM, RUMIDAS_GMX)
      tau = exp(m + theta(1) * c1)
    case (RUMIDAS_GM2M)
      tau = exp(m + theta(1) * c1 + theta(2) * c2)
    case (RUMIDAS_DAGM, RUMIDAS_DAGMX)
      tau = exp(m + theta(1) * c1 + theta(2) * c2)
    case (RUMIDAS_DAGM2M)
      tau = exp(m + theta(1) * c1 + theta(2) * c2 + theta(3) * c3 + theta(4) * c4)
    end select
    if (any(.not. (tau > tiny(1.0_dp))) .or. any(tau >= huge(1.0_dp) / 100.0_dp)) then
      status = RUMIDAS_NUMERICAL_ERROR
      call allocate_empty_outputs(n, loglik, conditional_volatility, long_run_volatility, short_run_variance)
      return
    end if

    allocate(loglik(n), conditional_volatility(n), long_run_volatility(n), short_run_variance(n))
    short_run_variance = 1.0_dp
    intercept = 1.0_dp - alpha - beta - 0.5_dp * gamma_value
    do i = 1, n
      step(i) = intercept + (alpha + merge(gamma_value, 0.0_dp, daily_ret(i) < 0.0_dp)) * &
        daily_ret(i) ** 2 / tau(i)
      if (has_x) step(i) = step(i) + x_coefficient * x_variable(i)
    end do
    do i = 2, n
      short_run_variance(i) = step(i - 1) + beta * short_run_variance(i - 1)
      if (short_run_variance(i) <= tiny(1.0_dp)) then
        status = RUMIDAS_INVALID_PARAMETER
        loglik = -huge(1.0_dp)
        conditional_volatility = huge(1.0_dp)
        long_run_volatility = sqrt(tau)
        return
      end if
    end do

    long_run_volatility = sqrt(tau)
    conditional_volatility = sqrt(short_run_variance * tau)
    mean_return = sample_mean(daily_ret)
    do i = 1, n
      variance = short_run_variance(i) * tau(i)
      select case (spec%distribution)
      case (RUMIDAS_NORMAL)
        residual = daily_ret(i) - mean_return
        loglik(i) = -0.5_dp * (LOG_TWO_PI + log(variance) + residual * residual / variance)
      case (RUMIDAS_STUDENT_T)
        ! The original rumidas code omits the constant -0.5*log(pi).  It is
        ! retained here for likelihood compatibility; parameter estimates are
        ! unaffected because the omitted term is constant.
        residual = daily_ret(i)
        loglik(i) = log_gamma(0.5_dp * (degrees + 1.0_dp)) - log_gamma(0.5_dp * degrees) - &
          0.5_dp * log((degrees - 2.0_dp) * variance) - 0.5_dp * (degrees + 1.0_dp) * &
          log(1.0_dp + residual * residual / ((degrees - 2.0_dp) * variance))
      case default
        status = RUMIDAS_INVALID_INPUT
        loglik = -huge(1.0_dp)
        return
      end select
    end do
  end subroutine garch_midas_evaluate

  pure logical function valid_garch_parameters(alpha, beta, gamma_value, skew, degrees, distribution)
    real(dp), intent(in) :: alpha, beta, gamma_value, degrees
    logical, intent(in) :: skew
    integer, intent(in) :: distribution
    valid_garch_parameters = alpha > 0.0_dp .and. beta > 0.0_dp
    if (skew) then
      valid_garch_parameters = valid_garch_parameters .and. alpha + beta + 0.5_dp * gamma_value < 1.0_dp
    else
      valid_garch_parameters = valid_garch_parameters .and. alpha + beta < 1.0_dp
    end if
    if (distribution == RUMIDAS_STUDENT_T) valid_garch_parameters = valid_garch_parameters .and. degrees > 2.0_dp
  end function valid_garch_parameters

  pure logical function valid_lag_parameters(w2, spec)
    real(dp), intent(in) :: w2(4)
    type(garch_midas_spec), intent(in) :: spec
    integer :: number, j
    select case (spec%model)
    case (RUMIDAS_GM, RUMIDAS_GMX)
      number = 1
    case (RUMIDAS_GM2M, RUMIDAS_DAGM, RUMIDAS_DAGMX)
      number = 2
    case (RUMIDAS_DAGM2M)
      number = 4
    case default
      valid_lag_parameters = .false.
      return
    end select
    valid_lag_parameters = .true.
    do j = 1, number
      if (spec%lag_function == RUMIDAS_BETA_LAG) valid_lag_parameters = valid_lag_parameters .and. w2(j) > 1.0_dp
      if (spec%lag_function == RUMIDAS_ALMON_LAG) valid_lag_parameters = valid_lag_parameters .and. w2(j) < 0.0_dp
    end do
  end function valid_lag_parameters

  subroutine allocate_empty_outputs(n, loglik, conditional, long_run, short_run)
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: loglik(:), conditional(:), long_run(:), short_run(:)
    allocate(loglik(max(n, 0)), conditional(max(n, 0)), long_run(max(n, 0)), short_run(max(n, 0)))
    loglik = -huge(1.0_dp)
    conditional = huge(1.0_dp)
    long_run = huge(1.0_dp)
    short_run = huge(1.0_dp)
  end subroutine allocate_empty_outputs

  subroutine evaluate_wrapper(param, spec, daily_ret, mv1, ll, cv, lr, status, mv2, x)
    real(dp), intent(in) :: param(:), daily_ret(:), mv1(:, :)
    type(garch_midas_spec), intent(in) :: spec
    real(dp), allocatable, intent(out), optional :: ll(:), cv(:), lr(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: mv2(:, :), x(:)
    real(dp), allocatable :: a(:), b(:), c(:), d(:)
    if (present(mv2) .and. present(x)) then
      call garch_midas_evaluate(param, spec, daily_ret, mv1, a, b, c, d, status, mv2, x)
    else if (present(mv2)) then
      call garch_midas_evaluate(param, spec, daily_ret, mv1, a, b, c, d, status, mv_m_2=mv2)
    else if (present(x)) then
      call garch_midas_evaluate(param, spec, daily_ret, mv1, a, b, c, d, status, x_variable=x)
    else
      call garch_midas_evaluate(param, spec, daily_ret, mv1, a, b, c, d, status)
    end if
    if (present(ll)) ll = a
    if (present(cv)) cv = b
    if (present(lr)) lr = c
  end subroutine evaluate_wrapper

  pure function make_spec(model, skew, k1, distribution, lag_function, k2) result(spec)
    integer, intent(in) :: model, k1
    logical, intent(in) :: skew
    integer, intent(in), optional :: distribution, lag_function, k2
    type(garch_midas_spec) :: spec
    spec%model = model
    spec%skew = skew
    spec%k1 = k1
    if (present(distribution)) spec%distribution = distribution
    if (present(lag_function)) spec%lag_function = lag_function
    if (present(k2)) spec%k2 = k2
  end function make_spec

  subroutine gm_loglik(param, daily_ret, mv_m, k, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM, .true., k, distribution, lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, ll=loglik, status=status)
  end subroutine gm_loglik

  subroutine gm_cond_vol(param, daily_ret, mv_m, k, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM, .true., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, cv=conditional, status=status)
  end subroutine gm_cond_vol

  subroutine gm_long_run_vol(param, daily_ret, mv_m, k, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM, .true., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, lr=long_run, status=status)
  end subroutine gm_long_run_vol


  subroutine gm_loglik_no_skew(param, daily_ret, mv_m, k, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM, .false., k, distribution, lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, ll=loglik, status=status)
  end subroutine gm_loglik_no_skew

  subroutine gm_cond_vol_no_skew(param, daily_ret, mv_m, k, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM, .false., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, cv=conditional, status=status)
  end subroutine gm_cond_vol_no_skew

  subroutine gm_long_run_vol_no_skew(param, daily_ret, mv_m, k, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM, .false., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, lr=long_run, status=status)
  end subroutine gm_long_run_vol_no_skew

  subroutine gm_x_loglik(param, daily_ret, x, mv_m, k, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GMX, .true., k, distribution, lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, ll=loglik, status=status, x=x)
  end subroutine gm_x_loglik

  subroutine gm_x_cond_vol(param, daily_ret, x, mv_m, k, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GMX, .true., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, cv=conditional, status=status, x=x)
  end subroutine gm_x_cond_vol

  subroutine gm_x_long_run_vol(param, daily_ret, x, mv_m, k, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GMX, .true., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, lr=long_run, status=status, x=x)
  end subroutine gm_x_long_run_vol

  subroutine gm_x_loglik_no_skew(param, daily_ret, x, mv_m, k, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GMX, .false., k, distribution, lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, ll=loglik, status=status, x=x)
  end subroutine gm_x_loglik_no_skew

  subroutine gm_x_cond_vol_no_skew(param, daily_ret, x, mv_m, k, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GMX, .false., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, cv=conditional, status=status, x=x)
  end subroutine gm_x_cond_vol_no_skew

  subroutine gm_x_long_run_vol_no_skew(param, daily_ret, x, mv_m, k, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GMX, .false., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, lr=long_run, status=status, x=x)
  end subroutine gm_x_long_run_vol_no_skew

  subroutine gm_2m_loglik(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM2M, .true., k_1, distribution, lag_function, k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, ll=loglik, status=status, mv2=mv_m_2)
  end subroutine gm_2m_loglik

  subroutine gm_2m_cond_vol(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM2M, .true., k_1, lag_function=lag_function, k2=k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, cv=conditional, status=status, mv2=mv_m_2)
  end subroutine gm_2m_cond_vol

  subroutine gm_2m_long_run_vol(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM2M, .true., k_1, lag_function=lag_function, k2=k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, lr=long_run, status=status, mv2=mv_m_2)
  end subroutine gm_2m_long_run_vol

  subroutine gm_2m_loglik_no_skew(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM2M, .false., k_1, distribution, lag_function, k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, ll=loglik, status=status, mv2=mv_m_2)
  end subroutine gm_2m_loglik_no_skew

  subroutine gm_2m_cond_vol_no_skew(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM2M, .false., k_1, lag_function=lag_function, k2=k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, cv=conditional, status=status, mv2=mv_m_2)
  end subroutine gm_2m_cond_vol_no_skew

  subroutine gm_2m_long_run_vol_no_skew(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_GM2M, .false., k_1, lag_function=lag_function, k2=k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, lr=long_run, status=status, mv2=mv_m_2)
  end subroutine gm_2m_long_run_vol_no_skew

  subroutine dagm_loglik(param, daily_ret, mv_m, k, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM, .true., k, distribution, lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, ll=loglik, status=status)
  end subroutine dagm_loglik

  subroutine dagm_cond_vol(param, daily_ret, mv_m, k, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM, .true., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, cv=conditional, status=status)
  end subroutine dagm_cond_vol

  subroutine dagm_long_run_vol(param, daily_ret, mv_m, k, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM, .true., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, lr=long_run, status=status)
  end subroutine dagm_long_run_vol

  subroutine dagm_loglik_no_skew(param, daily_ret, mv_m, k, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM, .false., k, distribution, lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, ll=loglik, status=status)
  end subroutine dagm_loglik_no_skew

  subroutine dagm_cond_vol_no_skew(param, daily_ret, mv_m, k, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM, .false., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, cv=conditional, status=status)
  end subroutine dagm_cond_vol_no_skew

  subroutine dagm_long_run_vol_no_skew(param, daily_ret, mv_m, k, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM, .false., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, lr=long_run, status=status)
  end subroutine dagm_long_run_vol_no_skew

  subroutine dagm_x_loglik(param, daily_ret, x, mv_m, k, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGMX, .true., k, distribution, lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, ll=loglik, status=status, x=x)
  end subroutine dagm_x_loglik

  subroutine dagm_x_cond_vol(param, daily_ret, x, mv_m, k, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGMX, .true., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, cv=conditional, status=status, x=x)
  end subroutine dagm_x_cond_vol

  subroutine dagm_x_long_run_vol(param, daily_ret, x, mv_m, k, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGMX, .true., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, lr=long_run, status=status, x=x)
  end subroutine dagm_x_long_run_vol

  subroutine dagm_x_loglik_no_skew(param, daily_ret, x, mv_m, k, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGMX, .false., k, distribution, lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, ll=loglik, status=status, x=x)
  end subroutine dagm_x_loglik_no_skew

  subroutine dagm_x_cond_vol_no_skew(param, daily_ret, x, mv_m, k, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGMX, .false., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, cv=conditional, status=status, x=x)
  end subroutine dagm_x_cond_vol_no_skew

  subroutine dagm_x_long_run_vol_no_skew(param, daily_ret, x, mv_m, k, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGMX, .false., k, lag_function=lag_function)
    call evaluate_wrapper(param, spec, daily_ret, mv_m, lr=long_run, status=status, x=x)
  end subroutine dagm_x_long_run_vol_no_skew

  subroutine dagm_2m_loglik(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM2M, .true., k_1, distribution, lag_function, k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, ll=loglik, status=status, mv2=mv_m_2)
  end subroutine dagm_2m_loglik

  subroutine dagm_2m_cond_vol(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM2M, .true., k_1, lag_function=lag_function, k2=k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, cv=conditional, status=status, mv2=mv_m_2)
  end subroutine dagm_2m_cond_vol

  subroutine dagm_2m_long_run(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM2M, .true., k_1, lag_function=lag_function, k2=k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, lr=long_run, status=status, mv2=mv_m_2)
  end subroutine dagm_2m_long_run

  subroutine dagm_2m_loglik_no_skew(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, loglik, status, distribution, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: distribution, lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM2M, .false., k_1, distribution, lag_function, k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, ll=loglik, status=status, mv2=mv_m_2)
  end subroutine dagm_2m_loglik_no_skew

  subroutine dagm_2m_cond_vol_no_skew(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, conditional, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: conditional(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM2M, .false., k_1, lag_function=lag_function, k2=k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, cv=conditional, status=status, mv2=mv_m_2)
  end subroutine dagm_2m_cond_vol_no_skew

  subroutine dagm_2m_long_run_no_skew(param, daily_ret, mv_m_1, mv_m_2, k_1, k_2, long_run, status, lag_function)
    real(dp), intent(in) :: param(:), daily_ret(:), mv_m_1(:, :), mv_m_2(:, :)
    integer, intent(in) :: k_1, k_2
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: lag_function
    type(garch_midas_spec) :: spec
    spec = make_spec(RUMIDAS_DAGM2M, .false., k_1, lag_function=lag_function, k2=k_2)
    call evaluate_wrapper(param, spec, daily_ret, mv_m_1, lr=long_run, status=status, mv2=mv_m_2)
  end subroutine dagm_2m_long_run_no_skew

end module rumidas_garch_midas
