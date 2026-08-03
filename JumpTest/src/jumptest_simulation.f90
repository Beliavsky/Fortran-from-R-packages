! SPDX-License-Identifier: MIT
module jumptest_simulation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use jumptest_kinds, only : dp, i8
  use jumptest_rng, only : rng_state, rng_seed, rng_normal, rng_poisson, rng_chisq
  use jumptest_status, only : JT_SUCCESS, JT_INVALID_ARGUMENT, JT_INVALID_DIMENSION, &
    JT_NONFINITE_INPUT, JT_NUMERICAL_FAILURE
  implicit none
  private

  type, public :: simulation_result
    real(dp), allocatable :: price(:)
    real(dp), allocatable :: variance(:,:)
    real(dp), allocatable :: jump_size(:)
    real(dp), allocatable :: jump_measure(:)
    integer, allocatable :: jump_count(:)
    integer :: status = JT_SUCCESS
    integer(i8) :: seed = 0_i8
  end type simulation_result

  public :: sv, svj, sv1f, sv1fj, sv2f
  public :: lp_path, pvc_path, pv2_path

contains

  subroutine sv(intervals_per_period, periods, result, p0, mu, v0, b, alpha, sigma, seed)
    integer, intent(in) :: intervals_per_period, periods
    type(simulation_result), intent(out) :: result
    real(dp), intent(in), optional :: p0, mu, v0, b, alpha, sigma
    integer(i8), intent(in), optional :: seed
    real(dp) :: p0v, muv, v0v, bv, alphav, sigmav, dt, df, c0
    real(dp), allocatable :: zvol(:), zprice(:), xchi(:), jump(:)
    real(dp), allocatable :: full_price(:), full_variance(:)
    type(rng_state) :: rng
    integer :: n, i

    call resolve_cir_arguments(intervals_per_period, periods, p0, mu, v0, b, alpha, &
      sigma, p0v, muv, v0v, bv, alphav, sigmav, n, result%status)
    if (result%status /= JT_SUCCESS) return
    call initialize_rng(rng, result, seed)
    dt = 1.0_dp/real(intervals_per_period, dp)
    df = 4.0_dp*bv*alphav/(sigmav*sigmav) - 1.0_dp
    if (df <= 0.0_dp) then
      result%status = JT_INVALID_ARGUMENT
      return
    end if
    c0 = cir_scale(alphav, sigmav, dt)
    allocate(zvol(n), zprice(n), xchi(n), jump(n))
    do i = 1, n
      zvol(i) = rng_normal(rng)
      zprice(i) = rng_normal(rng)
      xchi(i) = rng_chisq(rng, df)
    end do
    jump = 0.0_dp
    call lp_path(p0v, muv, v0v, dt, alphav, c0, zvol, zprice, jump, xchi, &
      full_price, full_variance, result%status)
    if (result%status /= JT_SUCCESS) return
    call move_alloc(full_price, result%price)
    allocate(result%variance(size(full_variance), 1))
    result%variance(:, 1) = full_variance
  end subroutine sv

  subroutine svj(intervals_per_period, periods, result, p0, lambda, mu, v0, b, &
      alpha, sigma, sigma1, seed)
    integer, intent(in) :: intervals_per_period, periods
    type(simulation_result), intent(out) :: result
    real(dp), intent(in), optional :: p0, lambda, mu, v0, b, alpha, sigma, sigma1
    integer(i8), intent(in), optional :: seed
    real(dp) :: p0v, lambdav, muv, v0v, bv, alphav, sigmav, sigma1v, dt, df, c0
    real(dp), allocatable :: zvol(:), zprice(:), xchi(:), jump(:)
    real(dp), allocatable :: full_price(:), full_variance(:)
    type(rng_state) :: rng
    integer :: n, i

    call resolve_cir_arguments(intervals_per_period, periods, p0, mu, v0, b, alpha, &
      sigma, p0v, muv, v0v, bv, alphav, sigmav, n, result%status)
    if (result%status /= JT_SUCCESS) return
    lambdav = 0.2_dp
    sigma1v = 1.0_dp
    if (present(lambda)) lambdav = lambda
    if (present(sigma1)) sigma1v = sigma1
    if (.not. ieee_is_finite(lambdav) .or. .not. ieee_is_finite(sigma1v) .or. &
        lambdav < 0.0_dp .or. sigma1v < 0.0_dp) then
      result%status = JT_INVALID_ARGUMENT
      return
    end if
    call initialize_rng(rng, result, seed)
    dt = 1.0_dp/real(intervals_per_period, dp)
    df = 4.0_dp*bv*alphav/(sigmav*sigmav) - 1.0_dp
    if (df <= 0.0_dp) then
      result%status = JT_INVALID_ARGUMENT
      return
    end if
    c0 = cir_scale(alphav, sigmav, dt)
    allocate(zvol(n), zprice(n), xchi(n), jump(n), result%jump_count(n), &
      result%jump_measure(n), result%jump_size(n))
    do i = 1, n
      zvol(i) = rng_normal(rng)
      zprice(i) = rng_normal(rng)
      xchi(i) = rng_chisq(rng, df)
      result%jump_count(i) = rng_poisson(rng, lambdav*dt)
      result%jump_measure(i) = sqrt(real(result%jump_count(i), dp)*sigma1v)
      jump(i) = result%jump_measure(i)*rng_normal(rng)
    end do
    result%jump_size = jump
    call lp_path(p0v, muv, v0v, dt, alphav, c0, zvol, zprice, jump, xchi, &
      full_price, full_variance, result%status)
    if (result%status /= JT_SUCCESS) return
    allocate(result%price(n), result%variance(n, 1))
    result%price = full_price(2:n + 1)
    result%variance(:, 1) = full_variance(2:n + 1)
  end subroutine svj

  subroutine sv1f(intervals_per_period, periods, result, p0, mu, v0, beta0, beta1, &
      alphav, correlation, seed)
    integer, intent(in) :: intervals_per_period, periods
    type(simulation_result), intent(out) :: result
    real(dp), intent(in), optional :: p0, mu, v0, beta0, beta1, alphav, correlation
    integer(i8), intent(in), optional :: seed
    real(dp) :: p0v, muv, v0v, beta0v, beta1v, alphavv, corrv, dt
    real(dp), allocatable :: shocks(:,:), jump(:), full_price(:), full_variance(:)
    type(rng_state) :: rng
    integer :: n, i

    call resolve_one_factor_arguments(intervals_per_period, periods, p0, mu, v0, beta0, &
      beta1, alphav, correlation, p0v, muv, v0v, beta0v, beta1v, alphavv, corrv, &
      n, result%status)
    if (result%status /= JT_SUCCESS) return
    call initialize_rng(rng, result, seed)
    dt = 1.0_dp/real(intervals_per_period, dp)
    allocate(shocks(n, 2), jump(n))
    do i = 1, n
      call correlated_pair(rng, corrv, shocks(i, 1), shocks(i, 2))
    end do
    jump = 0.0_dp
    call pvc_path(p0v, muv*dt, beta0v, beta1v, v0v, sqrt(dt), &
      1.0_dp + alphavv*dt, shocks, jump, full_price, full_variance, result%status)
    if (result%status /= JT_SUCCESS) return
    allocate(result%price(n), result%variance(n, 1))
    result%price = full_price(2:n + 1)
    result%variance(:, 1) = full_variance(2:n + 1)
  end subroutine sv1f

  subroutine sv1fj(intervals_per_period, periods, result, p0, lambda, mu, v0, beta0, &
      beta1, alphav, correlation, seed)
    integer, intent(in) :: intervals_per_period, periods
    type(simulation_result), intent(out) :: result
    real(dp), intent(in), optional :: p0, lambda, mu, v0, beta0, beta1, alphav, correlation
    integer(i8), intent(in), optional :: seed
    real(dp) :: p0v, lambdav, muv, v0v, beta0v, beta1v, alphavv, corrv, dt
    real(dp), allocatable :: shocks(:,:), jump(:), full_price(:), full_variance(:)
    type(rng_state) :: rng
    integer :: n, i

    call resolve_one_factor_arguments(intervals_per_period, periods, p0, mu, v0, beta0, &
      beta1, alphav, correlation, p0v, muv, v0v, beta0v, beta1v, alphavv, corrv, &
      n, result%status)
    if (result%status /= JT_SUCCESS) return
    lambdav = 0.2_dp
    if (present(lambda)) lambdav = lambda
    if (.not. ieee_is_finite(lambdav) .or. lambdav < 0.0_dp) then
      result%status = JT_INVALID_ARGUMENT
      return
    end if
    call initialize_rng(rng, result, seed)
    dt = 1.0_dp/real(intervals_per_period, dp)
    allocate(shocks(n, 2), jump(n), result%jump_count(n), result%jump_measure(n), &
      result%jump_size(n))
    do i = 1, n
      call correlated_pair(rng, corrv, shocks(i, 1), shocks(i, 2))
      result%jump_count(i) = rng_poisson(rng, lambdav*dt)
      result%jump_measure(i) = real(result%jump_count(i), dp)
      jump(i) = result%jump_measure(i)*rng_normal(rng)
    end do
    result%jump_size = jump
    call pvc_path(p0v, muv*dt, beta0v, beta1v, v0v, sqrt(dt), &
      1.0_dp + alphavv*dt, shocks, jump, full_price, full_variance, result%status)
    if (result%status /= JT_SUCCESS) return
    allocate(result%price(n), result%variance(n, 1))
    result%price = full_price(2:n + 1)
    result%variance(:, 1) = full_variance(2:n + 1)
  end subroutine sv1fj

  subroutine sv2f(intervals_per_period, periods, result, p0, mu, v1, v2, beta0, &
      beta1, beta2, alpha1, alpha2, beta_v2, r1, r2, seed)
    integer, intent(in) :: intervals_per_period, periods
    type(simulation_result), intent(out) :: result
    real(dp), intent(in), optional :: p0, mu, v1, v2, beta0, beta1, beta2
    real(dp), intent(in), optional :: alpha1, alpha2, beta_v2, r1, r2
    integer(i8), intent(in), optional :: seed
    real(dp) :: p0v, muv, v1v, v2v, beta0v, beta1v, beta2v
    real(dp) :: alpha1v, alpha2v, beta_v2v, r1v, r2v, dt
    real(dp) :: covariance(3, 3), lower(3, 3), independent(3)
    real(dp), allocatable :: shocks(:,:), full_price(:), full_variance(:,:)
    type(rng_state) :: rng
    integer :: n, i

    call resolve_two_factor_arguments(intervals_per_period, periods, p0, mu, v1, v2, &
      beta0, beta1, beta2, alpha1, alpha2, beta_v2, r1, r2, p0v, muv, v1v, v2v, &
      beta0v, beta1v, beta2v, alpha1v, alpha2v, beta_v2v, r1v, r2v, n, result%status)
    if (result%status /= JT_SUCCESS) return
    covariance = reshape([1.0_dp, r1v, r2v, r1v, 1.0_dp, 0.0_dp, &
      r2v, 0.0_dp, 1.0_dp], [3, 3])
    call cholesky_semidefinite(covariance, lower, result%status)
    if (result%status /= JT_SUCCESS) return
    call initialize_rng(rng, result, seed)
    allocate(shocks(n, 3))
    do i = 1, n
      independent = [rng_normal(rng), rng_normal(rng), rng_normal(rng)]
      shocks(i, :) = matmul(lower, independent)
    end do
    dt = 1.0_dp/real(intervals_per_period, dp)
    call pv2_path(p0v, muv*dt, beta0v, beta1v, beta2v, v1v, v2v, sqrt(dt), &
      1.0_dp + alpha1v*dt, 1.0_dp + alpha2v*dt, beta_v2v, shocks, &
      full_price, full_variance, result%status)
    if (result%status /= JT_SUCCESS) return
    allocate(result%price(n), result%variance(n, 2))
    result%price = full_price(2:n + 1)
    result%variance = full_variance(2:n + 1, :)
  end subroutine sv2f

  subroutine lp_path(p0, mu, v0, dt, alpha, c0, zvol, zprice, jump, xchi, price, &
      variance, status)
    real(dp), intent(in) :: p0, mu, v0, dt, alpha, c0
    real(dp), intent(in) :: zvol(:), zprice(:), jump(:), xchi(:)
    real(dp), allocatable, intent(out) :: price(:), variance(:)
    integer, intent(out) :: status
    integer :: i, n
    real(dp) :: variance_old

    n = size(zvol)
    if (size(zprice) /= n .or. size(jump) /= n .or. size(xchi) /= n) then
      status = JT_INVALID_DIMENSION
      allocate(price(0), variance(0))
      return
    end if
    if (dt <= 0.0_dp .or. c0 <= 0.0_dp .or. v0 < 0.0_dp .or. &
        .not. all(ieee_is_finite([p0, mu, v0, dt, alpha, c0])) .or. &
        .not. all(ieee_is_finite(zvol)) .or. .not. all(ieee_is_finite(zprice)) .or. &
        .not. all(ieee_is_finite(jump)) .or. .not. all(ieee_is_finite(xchi))) then
      status = JT_INVALID_ARGUMENT
      allocate(price(0), variance(0))
      return
    end if
    allocate(price(n + 1), variance(n + 1))
    price(1) = p0
    variance(1) = v0
    do i = 1, n
      variance_old = max(0.0_dp, variance(i))
      price(i + 1) = price(i) + (mu - 0.5_dp*variance_old)*dt + &
        sqrt(variance_old*dt)*zprice(i) + jump(i)
      variance(i + 1) = c0*((zvol(i) + sqrt(variance_old/c0)* &
        exp(-0.5_dp*alpha*dt))**2 + max(0.0_dp, xchi(i)))
      if (.not. ieee_is_finite(price(i + 1)) .or. .not. ieee_is_finite(variance(i + 1))) then
        status = JT_NUMERICAL_FAILURE
        return
      end if
    end do
    status = JT_SUCCESS
  end subroutine lp_path

  subroutine pvc_path(p0, mean_step, beta0, beta1, v0, shock_scale, persistence, &
      shocks, jump, price, variance, status)
    real(dp), intent(in) :: p0, mean_step, beta0, beta1, v0, shock_scale, persistence
    real(dp), intent(in) :: shocks(:,:), jump(:)
    real(dp), allocatable, intent(out) :: price(:), variance(:)
    integer, intent(out) :: status
    integer :: i, n

    n = size(shocks, 1)
    if (size(shocks, 2) /= 2 .or. size(jump) /= n) then
      status = JT_INVALID_DIMENSION
      allocate(price(0), variance(0))
      return
    end if
    if (shock_scale < 0.0_dp .or. .not. all(ieee_is_finite(shocks)) .or. &
        .not. all(ieee_is_finite(jump)) .or. .not. all(ieee_is_finite( &
        [p0, mean_step, beta0, beta1, v0, shock_scale, persistence]))) then
      status = JT_INVALID_ARGUMENT
      allocate(price(0), variance(0))
      return
    end if
    allocate(price(n + 1), variance(n + 1))
    price(1) = p0
    variance(1) = v0
    do i = 1, n
      price(i + 1) = price(i) + mean_step + exp(beta0 + beta1*variance(i))* &
        shock_scale*shocks(i, 1) + jump(i)
      variance(i + 1) = persistence*variance(i) + shock_scale*shocks(i, 2)
      if (.not. ieee_is_finite(price(i + 1)) .or. .not. ieee_is_finite(variance(i + 1))) then
        status = JT_NUMERICAL_FAILURE
        return
      end if
    end do
    status = JT_SUCCESS
  end subroutine pvc_path

  subroutine pv2_path(p0, mean_step, beta0, beta1, beta2, v10, v20, shock_scale, &
      persistence1, persistence2, beta_v2, shocks, price, variance, status)
    real(dp), intent(in) :: p0, mean_step, beta0, beta1, beta2, v10, v20
    real(dp), intent(in) :: shock_scale, persistence1, persistence2, beta_v2
    real(dp), intent(in) :: shocks(:,:)
    real(dp), allocatable, intent(out) :: price(:), variance(:,:)
    integer, intent(out) :: status
    integer :: i, n

    n = size(shocks, 1)
    if (size(shocks, 2) /= 3) then
      status = JT_INVALID_DIMENSION
      allocate(price(0), variance(0, 0))
      return
    end if
    if (shock_scale < 0.0_dp .or. .not. all(ieee_is_finite(shocks)) .or. &
        .not. all(ieee_is_finite([p0, mean_step, beta0, beta1, beta2, v10, v20, &
        shock_scale, persistence1, persistence2, beta_v2]))) then
      status = JT_INVALID_ARGUMENT
      allocate(price(0), variance(0, 0))
      return
    end if
    allocate(price(n + 1), variance(n + 1, 2))
    price(1) = p0
    variance(1, :) = [v10, v20]
    do i = 1, n
      price(i + 1) = price(i) + mean_step - exp(beta0 + beta1*variance(i, 1) + &
        beta2*variance(i, 2))*shock_scale*shocks(i, 1)
      variance(i + 1, 1) = persistence1*variance(i, 1) + shock_scale*shocks(i, 2)
      variance(i + 1, 2) = persistence2*variance(i, 2) + &
        (1.0_dp + beta_v2*variance(i, 2))*shock_scale*shocks(i, 3)
      if (.not. ieee_is_finite(price(i + 1)) .or. &
          .not. all(ieee_is_finite(variance(i + 1, :)))) then
        status = JT_NUMERICAL_FAILURE
        return
      end if
    end do
    status = JT_SUCCESS
  end subroutine pv2_path

  subroutine initialize_rng(rng, result, seed)
    type(rng_state), intent(out) :: rng
    type(simulation_result), intent(inout) :: result
    integer(i8), intent(in), optional :: seed

    result%seed = 104729_i8
    if (present(seed)) result%seed = seed
    call rng_seed(rng, result%seed)
  end subroutine initialize_rng

  subroutine resolve_cir_arguments(intervals, periods, p0, mu, v0, b, alpha, sigma, &
      p0v, muv, v0v, bv, alphav, sigmav, n, status)
    integer, intent(in) :: intervals, periods
    real(dp), intent(in), optional :: p0, mu, v0, b, alpha, sigma
    real(dp), intent(out) :: p0v, muv, v0v, bv, alphav, sigmav
    integer, intent(out) :: n, status
    integer(i8) :: n64

    p0v = 3.0_dp; muv = 0.05_dp; v0v = 0.0_dp
    bv = 0.2_dp; alphav = 0.015_dp; sigmav = 0.05_dp
    if (present(p0)) p0v = p0
    if (present(mu)) muv = mu
    if (present(v0)) v0v = v0
    if (present(b)) bv = b
    if (present(alpha)) alphav = alpha
    if (present(sigma)) sigmav = sigma
    call validate_count(intervals, periods, n64, status)
    if (status /= JT_SUCCESS) then
      n = 0
      return
    end if
    n = int(n64)
    if (.not. all(ieee_is_finite([p0v, muv, v0v, bv, alphav, sigmav])) .or. &
        v0v < 0.0_dp .or. bv <= 0.0_dp .or. alphav <= 0.0_dp .or. sigmav <= 0.0_dp) then
      status = JT_INVALID_ARGUMENT
    end if
  end subroutine resolve_cir_arguments

  subroutine resolve_one_factor_arguments(intervals, periods, p0, mu, v0, beta0, beta1, &
      alphav, correlation, p0v, muv, v0v, beta0v, beta1v, alphavv, corrv, n, status)
    integer, intent(in) :: intervals, periods
    real(dp), intent(in), optional :: p0, mu, v0, beta0, beta1, alphav, correlation
    real(dp), intent(out) :: p0v, muv, v0v, beta0v, beta1v, alphavv, corrv
    integer, intent(out) :: n, status
    integer(i8) :: n64

    p0v = 3.0_dp; muv = 0.03_dp; v0v = 5.0_dp
    beta0v = 0.0_dp; beta1v = 0.125_dp; alphavv = -0.1_dp; corrv = -0.62_dp
    if (present(p0)) p0v = p0
    if (present(mu)) muv = mu
    if (present(v0)) v0v = v0
    if (present(beta0)) beta0v = beta0
    if (present(beta1)) beta1v = beta1
    if (present(alphav)) alphavv = alphav
    if (present(correlation)) corrv = correlation
    call validate_count(intervals, periods, n64, status)
    if (status /= JT_SUCCESS) then
      n = 0
      return
    end if
    n = int(n64)
    if (.not. all(ieee_is_finite([p0v, muv, v0v, beta0v, beta1v, alphavv, corrv])) .or. &
        abs(corrv) > 1.0_dp) status = JT_INVALID_ARGUMENT
  end subroutine resolve_one_factor_arguments

  subroutine resolve_two_factor_arguments(intervals, periods, p0, mu, v1, v2, beta0, &
      beta1, beta2, alpha1, alpha2, beta_v2, r1, r2, p0v, muv, v1v, v2v, beta0v, &
      beta1v, beta2v, alpha1v, alpha2v, beta_v2v, r1v, r2v, n, status)
    integer, intent(in) :: intervals, periods
    real(dp), intent(in), optional :: p0, mu, v1, v2, beta0, beta1, beta2
    real(dp), intent(in), optional :: alpha1, alpha2, beta_v2, r1, r2
    real(dp), intent(out) :: p0v, muv, v1v, v2v, beta0v, beta1v, beta2v
    real(dp), intent(out) :: alpha1v, alpha2v, beta_v2v, r1v, r2v
    integer, intent(out) :: n, status
    integer(i8) :: n64

    p0v = 3.0_dp; muv = 0.03_dp; v1v = 0.5_dp; v2v = 0.5_dp
    beta0v = -1.2_dp; beta1v = 0.04_dp; beta2v = 1.5_dp
    alpha1v = -0.137_dp*exp(-2.0_dp); alpha2v = -1.386_dp
    beta_v2v = 0.25_dp; r1v = -0.3_dp; r2v = -0.3_dp
    if (present(p0)) p0v = p0
    if (present(mu)) muv = mu
    if (present(v1)) v1v = v1
    if (present(v2)) v2v = v2
    if (present(beta0)) beta0v = beta0
    if (present(beta1)) beta1v = beta1
    if (present(beta2)) beta2v = beta2
    if (present(alpha1)) alpha1v = alpha1
    if (present(alpha2)) alpha2v = alpha2
    if (present(beta_v2)) beta_v2v = beta_v2
    if (present(r1)) r1v = r1
    if (present(r2)) r2v = r2
    call validate_count(intervals, periods, n64, status)
    if (status /= JT_SUCCESS) then
      n = 0
      return
    end if
    n = int(n64)
    if (.not. all(ieee_is_finite([p0v, muv, v1v, v2v, beta0v, beta1v, beta2v, &
        alpha1v, alpha2v, beta_v2v, r1v, r2v])) .or. abs(r1v) > 1.0_dp .or. &
        abs(r2v) > 1.0_dp) status = JT_INVALID_ARGUMENT
  end subroutine resolve_two_factor_arguments

  subroutine validate_count(intervals, periods, n64, status)
    integer, intent(in) :: intervals, periods
    integer(i8), intent(out) :: n64
    integer, intent(out) :: status

    if (intervals < 1 .or. periods < 1) then
      n64 = 0_i8
      status = JT_INVALID_DIMENSION
      return
    end if
    n64 = int(intervals, i8)*int(periods, i8)
    if (n64 > int(huge(0), i8)) then
      status = JT_INVALID_DIMENSION
    else
      status = JT_SUCCESS
    end if
  end subroutine validate_count

  pure function cir_scale(alpha, sigma, dt) result(c0)
    real(dp), intent(in) :: alpha, sigma, dt
    real(dp) :: c0

    if (abs(alpha*dt) < 1.0e-8_dp) then
      c0 = 0.25_dp*sigma*sigma*dt
    else
      c0 = 0.25_dp*sigma*sigma*(1.0_dp - exp(-alpha*dt))/alpha
    end if
  end function cir_scale

  subroutine correlated_pair(rng, correlation, first, second)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: correlation
    real(dp), intent(out) :: first, second
    real(dp) :: independent

    first = rng_normal(rng)
    independent = rng_normal(rng)
    second = correlation*first + sqrt(max(0.0_dp, 1.0_dp - correlation*correlation))*independent
  end subroutine correlated_pair

  subroutine cholesky_semidefinite(matrix, lower, status)
    real(dp), intent(in) :: matrix(:,:)
    real(dp), intent(out) :: lower(:,:)
    integer, intent(out) :: status
    real(dp) :: value, tolerance
    integer :: i, j, k, n

    n = size(matrix, 1)
    if (size(matrix, 2) /= n .or. size(lower, 1) /= n .or. size(lower, 2) /= n) then
      status = JT_INVALID_DIMENSION
      return
    end if
    lower = 0.0_dp
    tolerance = 1.0e-12_dp*max(1.0_dp, maxval(abs(matrix)))
    do i = 1, n
      do j = 1, i
        value = matrix(i, j)
        do k = 1, j - 1
          value = value - lower(i, k)*lower(j, k)
        end do
        if (i == j) then
          if (value < -tolerance) then
            status = JT_INVALID_ARGUMENT
            return
          end if
          lower(i, j) = sqrt(max(0.0_dp, value))
        else if (lower(j, j) > tolerance) then
          lower(i, j) = value/lower(j, j)
        else if (abs(value) > tolerance) then
          status = JT_INVALID_ARGUMENT
          return
        end if
      end do
    end do
    status = JT_SUCCESS
  end subroutine cholesky_semidefinite

end module jumptest_simulation
