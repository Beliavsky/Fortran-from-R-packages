! SPDX-License-Identifier: GPL-3.0-only
module rumidas_mem_models
  use rumidas_kinds, only: dp
  use rumidas_status
  use rumidas_types
  use rumidas_weights, only: midas_weighted_component
  use rumidas_statistics, only: sample_mean
  implicit none
  private

  public :: mem_parameter_count, mem_evaluate
  public :: mem_loglik, mem_pred, mem_loglik_no_skew, mem_pred_no_skew
  public :: mem_x_loglik, mem_x_pred, mem_x_loglik_no_skew, mem_x_pred_no_skew
  public :: mem_midas_loglik, mem_midas_pred, mem_midas_lr_pred
  public :: mem_midas_loglik_no_skew, mem_midas_pred_no_skew, mem_midas_lr_pred_no_skew
  public :: mem_midas_x_loglik, mem_midas_x_pred, mem_midas_x_lr_pred
  public :: mem_midas_x_loglik_no_skew, mem_midas_x_pred_no_skew, mem_midas_x_lr_pred_no_skew

contains

  pure integer function mem_parameter_count(spec) result(number)
    type(mem_spec), intent(in) :: spec
    number = 2
    if (spec%skew) number = number + 1
    select case (spec%model)
    case (RUMIDAS_MEM)
    case (RUMIDAS_MEM_X)
      number = number + 1
    case (RUMIDAS_MEM_MIDAS)
      number = number + 3
    case (RUMIDAS_MEM_MIDAS_X)
      number = number + 4
    case default
      number = 0
    end select
  end function mem_parameter_count

  subroutine mem_evaluate(param, spec, x, loglik, prediction, long_run, short_run, status, &
      daily_ret, mv_m, z_variable)
    real(dp), intent(in) :: param(:), x(:)
    type(mem_spec), intent(in) :: spec
    real(dp), allocatable, intent(out) :: loglik(:), prediction(:), long_run(:), short_run(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: daily_ret(:), mv_m(:, :), z_variable(:)

    real(dp), allocatable :: tau(:), component(:), step(:)
    real(dp) :: alpha, beta, gamma_value, m, theta, w2, delta, mean_x, intercept
    integer :: n, i, index, local_status
    logical :: has_midas, has_x

    status = RUMIDAS_SUCCESS
    n = size(x)
    call allocate_mem_outputs(n, loglik, prediction, long_run, short_run)
    if (n <= 0 .or. size(param) /= mem_parameter_count(spec) .or. any(x <= 0.0_dp)) then
      status = RUMIDAS_INVALID_INPUT
      return
    end if
    has_midas = spec%model == RUMIDAS_MEM_MIDAS .or. spec%model == RUMIDAS_MEM_MIDAS_X
    has_x = spec%model == RUMIDAS_MEM_X .or. spec%model == RUMIDAS_MEM_MIDAS_X
    if (spec%skew) then
      if (.not. present(daily_ret)) then
        status = RUMIDAS_INVALID_INPUT
        return
      end if
      if (size(daily_ret) /= n) then
        status = RUMIDAS_DIMENSION_ERROR
        return
      end if
    end if
    if (has_midas) then
      if (.not. present(mv_m)) then
        status = RUMIDAS_INVALID_INPUT
        return
      end if
      if (size(mv_m, 1) /= spec%k + 1 .or. size(mv_m, 2) /= n) then
        status = RUMIDAS_DIMENSION_ERROR
        return
      end if
    end if
    if (has_x) then
      if (.not. present(z_variable)) then
        status = RUMIDAS_INVALID_INPUT
        return
      end if
      if (size(z_variable) /= n) then
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
    m = 0.0_dp; theta = 0.0_dp; w2 = 2.0_dp; delta = 0.0_dp
    if (has_midas) then
      m = param(index)
      theta = param(index + 1)
      w2 = param(index + 2)
      index = index + 3
    end if
    if (has_x) delta = param(index)

    if (alpha <= 0.0_dp .or. beta <= 0.0_dp .or. &
        alpha + beta + 0.5_dp * gamma_value >= 1.0_dp) then
      status = RUMIDAS_INVALID_PARAMETER
      return
    end if
    if (has_midas .and. w2 <= 1.0_dp) then
      status = RUMIDAS_INVALID_PARAMETER
      return
    end if

    allocate(tau(n), component(n), step(n))
    tau = 1.0_dp
    if (has_midas) then
      call midas_weighted_component(mv_m, spec%k, w2, RUMIDAS_BETA_LAG, component, local_status)
      if (local_status /= RUMIDAS_SUCCESS) then
        status = local_status
        return
      end if
      tau = exp(m + theta * component)
      if (any(tau <= tiny(1.0_dp)) .or. any(tau >= huge(1.0_dp) / 100.0_dp)) then
        status = RUMIDAS_NUMERICAL_ERROR
        return
      end if
    end if

    mean_x = sample_mean(x)
    short_run = mean_x
    intercept = (1.0_dp - alpha - beta - 0.5_dp * gamma_value) * mean_x
    do i = 1, n
      step(i) = intercept + alpha * x(i) / tau(i)
      if (spec%skew) then
        if (daily_ret(i) < 0.0_dp) step(i) = step(i) + gamma_value * x(i) / tau(i)
      end if
      if (has_x) step(i) = step(i) + delta * z_variable(i)
    end do
    do i = 2, n
      short_run(i) = step(i - 1) + beta * short_run(i - 1)
      if (short_run(i) <= tiny(1.0_dp)) then
        status = RUMIDAS_INVALID_PARAMETER
        return
      end if
    end do
    long_run = tau
    prediction = short_run * tau
    if (any(prediction <= tiny(1.0_dp))) then
      status = RUMIDAS_NUMERICAL_ERROR
      return
    end if
    loglik = -log(prediction) - x / prediction
  end subroutine mem_evaluate

  subroutine allocate_mem_outputs(n, loglik, prediction, long_run, short_run)
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: loglik(:), prediction(:), long_run(:), short_run(:)
    allocate(loglik(max(n, 0)), prediction(max(n, 0)), long_run(max(n, 0)), short_run(max(n, 0)))
    loglik = -huge(1.0_dp)
    prediction = huge(1.0_dp)
    long_run = huge(1.0_dp)
    short_run = huge(1.0_dp)
  end subroutine allocate_mem_outputs

  pure function make_mem_spec(model, skew, k) result(spec)
    integer, intent(in) :: model
    logical, intent(in) :: skew
    integer, intent(in), optional :: k
    type(mem_spec) :: spec
    spec%model = model
    spec%skew = skew
    if (present(k)) spec%k = k
  end function make_mem_spec

  subroutine mem_wrapper(param, spec, x, status, ll, pred, lr, daily_ret, mv_m, z)
    real(dp), intent(in) :: param(:), x(:)
    type(mem_spec), intent(in) :: spec
    integer, intent(out) :: status
    real(dp), allocatable, intent(out), optional :: ll(:), pred(:), lr(:)
    real(dp), intent(in), optional :: daily_ret(:), mv_m(:, :), z(:)
    real(dp), allocatable :: a(:), b(:), c(:), d(:)
    if (present(daily_ret) .and. present(mv_m) .and. present(z)) then
      call mem_evaluate(param, spec, x, a, b, c, d, status, daily_ret, mv_m, z)
    else if (present(daily_ret) .and. present(mv_m)) then
      call mem_evaluate(param, spec, x, a, b, c, d, status, daily_ret=daily_ret, mv_m=mv_m)
    else if (present(daily_ret) .and. present(z)) then
      call mem_evaluate(param, spec, x, a, b, c, d, status, daily_ret=daily_ret, z_variable=z)
    else if (present(mv_m) .and. present(z)) then
      call mem_evaluate(param, spec, x, a, b, c, d, status, mv_m=mv_m, z_variable=z)
    else if (present(daily_ret)) then
      call mem_evaluate(param, spec, x, a, b, c, d, status, daily_ret=daily_ret)
    else if (present(mv_m)) then
      call mem_evaluate(param, spec, x, a, b, c, d, status, mv_m=mv_m)
    else if (present(z)) then
      call mem_evaluate(param, spec, x, a, b, c, d, status, z_variable=z)
    else
      call mem_evaluate(param, spec, x, a, b, c, d, status)
    end if
    if (present(ll)) ll = a
    if (present(pred)) pred = b
    if (present(lr)) lr = c
  end subroutine mem_wrapper

  subroutine mem_loglik(param, x, daily_ret, loglik, status)
    real(dp), intent(in) :: param(:), x(:), daily_ret(:)
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM, .true.), x, status, ll=loglik, daily_ret=daily_ret)
  end subroutine mem_loglik

  subroutine mem_pred(param, x, daily_ret, prediction, status)
    real(dp), intent(in) :: param(:), x(:), daily_ret(:)
    real(dp), allocatable, intent(out) :: prediction(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM, .true.), x, status, pred=prediction, daily_ret=daily_ret)
  end subroutine mem_pred

  subroutine mem_loglik_no_skew(param, x, loglik, status)
    real(dp), intent(in) :: param(:), x(:)
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM, .false.), x, status, ll=loglik)
  end subroutine mem_loglik_no_skew

  subroutine mem_pred_no_skew(param, x, prediction, status)
    real(dp), intent(in) :: param(:), x(:)
    real(dp), allocatable, intent(out) :: prediction(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM, .false.), x, status, pred=prediction)
  end subroutine mem_pred_no_skew

  subroutine mem_x_loglik(param, x, daily_ret, z, loglik, status)
    real(dp), intent(in) :: param(:), x(:), daily_ret(:), z(:)
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_X, .true.), x, status, ll=loglik, daily_ret=daily_ret, z=z)
  end subroutine mem_x_loglik

  subroutine mem_x_pred(param, x, daily_ret, z, prediction, status)
    real(dp), intent(in) :: param(:), x(:), daily_ret(:), z(:)
    real(dp), allocatable, intent(out) :: prediction(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_X, .true.), x, status, pred=prediction, daily_ret=daily_ret, z=z)
  end subroutine mem_x_pred

  subroutine mem_x_loglik_no_skew(param, x, z, loglik, status)
    real(dp), intent(in) :: param(:), x(:), z(:)
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_X, .false.), x, status, ll=loglik, z=z)
  end subroutine mem_x_loglik_no_skew

  subroutine mem_x_pred_no_skew(param, x, z, prediction, status)
    real(dp), intent(in) :: param(:), x(:), z(:)
    real(dp), allocatable, intent(out) :: prediction(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_X, .false.), x, status, pred=prediction, z=z)
  end subroutine mem_x_pred_no_skew

  subroutine mem_midas_loglik(param, x, daily_ret, mv_m, k, loglik, status)
    real(dp), intent(in) :: param(:), x(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS, .true., k), x, status, ll=loglik, &
      daily_ret=daily_ret, mv_m=mv_m)
  end subroutine mem_midas_loglik

  subroutine mem_midas_pred(param, x, daily_ret, mv_m, k, prediction, status)
    real(dp), intent(in) :: param(:), x(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: prediction(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS, .true., k), x, status, pred=prediction, &
      daily_ret=daily_ret, mv_m=mv_m)
  end subroutine mem_midas_pred

  subroutine mem_midas_lr_pred(param, x, daily_ret, mv_m, k, long_run, status)
    real(dp), intent(in) :: param(:), x(:), daily_ret(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS, .true., k), x, status, lr=long_run, &
      daily_ret=daily_ret, mv_m=mv_m)
  end subroutine mem_midas_lr_pred

  subroutine mem_midas_loglik_no_skew(param, x, mv_m, k, loglik, status)
    real(dp), intent(in) :: param(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS, .false., k), x, status, ll=loglik, mv_m=mv_m)
  end subroutine mem_midas_loglik_no_skew

  subroutine mem_midas_pred_no_skew(param, x, mv_m, k, prediction, status)
    real(dp), intent(in) :: param(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: prediction(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS, .false., k), x, status, pred=prediction, mv_m=mv_m)
  end subroutine mem_midas_pred_no_skew

  subroutine mem_midas_lr_pred_no_skew(param, x, mv_m, k, long_run, status)
    real(dp), intent(in) :: param(:), x(:), mv_m(:, :)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS, .false., k), x, status, lr=long_run, mv_m=mv_m)
  end subroutine mem_midas_lr_pred_no_skew

  subroutine mem_midas_x_loglik(param, x, daily_ret, mv_m, k, z, loglik, status)
    real(dp), intent(in) :: param(:), x(:), daily_ret(:), mv_m(:, :), z(:)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS_X, .true., k), x, status, ll=loglik, &
      daily_ret=daily_ret, mv_m=mv_m, z=z)
  end subroutine mem_midas_x_loglik

  subroutine mem_midas_x_pred(param, x, daily_ret, mv_m, k, z, prediction, status)
    real(dp), intent(in) :: param(:), x(:), daily_ret(:), mv_m(:, :), z(:)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: prediction(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS_X, .true., k), x, status, pred=prediction, &
      daily_ret=daily_ret, mv_m=mv_m, z=z)
  end subroutine mem_midas_x_pred

  subroutine mem_midas_x_lr_pred(param, x, daily_ret, mv_m, k, z, long_run, status)
    real(dp), intent(in) :: param(:), x(:), daily_ret(:), mv_m(:, :), z(:)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS_X, .true., k), x, status, lr=long_run, &
      daily_ret=daily_ret, mv_m=mv_m, z=z)
  end subroutine mem_midas_x_lr_pred

  subroutine mem_midas_x_loglik_no_skew(param, x, mv_m, k, z, loglik, status)
    real(dp), intent(in) :: param(:), x(:), mv_m(:, :), z(:)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: loglik(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS_X, .false., k), x, status, ll=loglik, mv_m=mv_m, z=z)
  end subroutine mem_midas_x_loglik_no_skew

  subroutine mem_midas_x_pred_no_skew(param, x, mv_m, k, z, prediction, status)
    real(dp), intent(in) :: param(:), x(:), mv_m(:, :), z(:)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: prediction(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS_X, .false., k), x, status, pred=prediction, mv_m=mv_m, z=z)
  end subroutine mem_midas_x_pred_no_skew

  subroutine mem_midas_x_lr_pred_no_skew(param, x, mv_m, k, z, long_run, status)
    real(dp), intent(in) :: param(:), x(:), mv_m(:, :), z(:)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: long_run(:)
    integer, intent(out) :: status
    call mem_wrapper(param, make_mem_spec(RUMIDAS_MEM_MIDAS_X, .false., k), x, status, lr=long_run, mv_m=mv_m, z=z)
  end subroutine mem_midas_x_lr_pred_no_skew

end module rumidas_mem_models
