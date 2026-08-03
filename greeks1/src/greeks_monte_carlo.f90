! SPDX-License-Identifier: MIT
module greeks_monte_carlo
  use greeks_kinds, only: dp
  use greeks_types, only: greek_result, initialize_result, set_error
  use greeks_types, only: payoff_callback
  use greeks_types, only: greeks_invalid_argument, greeks_unknown_payoff
  use greeks_math, only: rng_state, rng_seed, rng_normal, rng_poisson
  use greeks_math, only: rng_student_t3, sample_mean, sample_stderr
  use greeks_math, only: regression_intercept
  use greeks_payoffs, only: payoff_value, payoff_derivative, valid_standard_payoff
  use greeks_integrals, only: calc_i, calc_i_1, calc_i_2, calc_i_3
  use greeks_integrals, only: calc_xw, calc_txw
  use greeks_black_scholes, only: bs_geometric_asian_greeks
  implicit none
  private
  public :: malliavin_european_greeks, malliavin_asian_greeks
  public :: malliavin_geometric_asian_greeks, bs_malliavin_asian_greeks

contains

  subroutine validate_mc(spot, strike, time, sigma, payoff, paths, steps, &
      result, valid, allow_custom)
    real(dp), intent(in) :: spot, strike, time, sigma
    character(len=*), intent(in) :: payoff
    integer, intent(in) :: paths, steps
    type(greek_result), intent(inout) :: result
    logical, intent(out) :: valid
    logical, intent(in), optional :: allow_custom
    logical :: custom
    custom = .false.
    if (present(allow_custom)) custom = allow_custom
    valid = .false.
    if (spot <= 0.0_dp .or. strike <= 0.0_dp .or. time <= 0.0_dp .or. &
        sigma <= 0.0_dp .or. paths < 2 .or. steps < 1) then
      call set_error(result, greeks_invalid_argument, &
        'invalid spot, strike, time, volatility, paths, or steps')
      return
    end if
    if (.not. valid_standard_payoff(payoff) .and. .not. custom) then
      call set_error(result, greeks_unknown_payoff, 'unsupported payoff')
      return
    end if
    valid = .true.
  end subroutine validate_mc

  subroutine simulate_paths(paths, steps, time, sigma, drift, seed, antithetic, &
      jump_diffusion, jump_intensity, jump_scale, w, x)
    integer, intent(in) :: paths, steps, seed
    real(dp), intent(in) :: time, sigma, drift, jump_intensity, jump_scale
    logical, intent(in) :: antithetic, jump_diffusion
    real(dp), allocatable, intent(out) :: w(:, :), x(:, :)
    type(rng_state) :: rng
    real(dp) :: dt, dz, jump_value
    integer :: i, j, pair_count, n_jumps, k

    dt = time/real(steps, dp)
    allocate(w(paths, steps + 1), x(paths, steps + 1))
    w(:, 1) = 0.0_dp
    x(:, 1) = 1.0_dp
    call rng_seed(rng, seed)
    if (antithetic) then
      pair_count = paths/2
      do j = 2, steps + 1
        do i = 1, pair_count
          dz = sqrt(dt)*rng_normal(rng)
          w(i, j) = w(i, j - 1) + dz
          w(i + pair_count, j) = w(i + pair_count, j - 1) - dz
        end do
        if (2*pair_count < paths) then
          dz = sqrt(dt)*rng_normal(rng)
          w(paths, j) = w(paths, j - 1) + dz
        end if
      end do
    else
      do j = 2, steps + 1
        do i = 1, paths
          w(i, j) = w(i, j - 1) + sqrt(dt)*rng_normal(rng)
        end do
      end do
    end if

    do j = 2, steps + 1
      x(:, j) = exp((drift - 0.5_dp*sigma**2)*real(j - 1, dp)*dt + &
        sigma*w(:, j))
    end do
    if (jump_diffusion) then
      do i = 1, paths
        jump_value = 0.0_dp
        do j = 2, steps + 1
          n_jumps = rng_poisson(rng, jump_intensity*dt)
          do k = 1, n_jumps
            jump_value = jump_value + jump_scale*rng_student_t3(rng)
          end do
          x(i, j) = x(i, j)*exp(jump_value)
        end do
      end do
    end if
  end subroutine simulate_paths

  subroutine weighted_estimate(samples, discount, result, index_value)
    real(dp), intent(in) :: samples(:), discount
    type(greek_result), intent(inout) :: result
    integer, intent(in) :: index_value
    result%values(index_value) = discount*sample_mean(samples)
    result%standard_errors(index_value) = discount*sample_stderr(samples)
  end subroutine weighted_estimate

  subroutine malliavin_european_greeks(spot, strike, rate, time, sigma, payoff, &
      requested, result, paths, seed, antithetic, payoff_fn)
    real(dp), intent(in) :: spot, strike, rate, time, sigma
    character(len=*), intent(in) :: payoff
    character(len=*), intent(in) :: requested(:)
    type(greek_result), intent(out) :: result
    integer, intent(in), optional :: paths, seed
    logical, intent(in), optional :: antithetic
    procedure(payoff_callback), optional :: payoff_fn
    integer :: npaths, seed_value, i, j
    logical :: anti, valid
    type(rng_state) :: rng
    real(dp), allocatable :: wt(:), xt(:), samples(:)
    real(dp) :: z

    call initialize_result(result, requested)
    npaths = 10000; if (present(paths)) npaths = paths
    seed_value = 1; if (present(seed)) seed_value = seed
    anti = .false.; if (present(antithetic)) anti = antithetic
    call validate_mc(spot, strike, time, sigma, payoff, npaths, 1, result, valid, &
      present(payoff_fn))
    if (.not. valid) return
    allocate(wt(npaths), xt(npaths), samples(npaths))
    call rng_seed(rng, seed_value)
    if (anti) then
      do i = 1, npaths/2
        z = sqrt(time)*rng_normal(rng)
        wt(i) = z
        wt(i + npaths/2) = -z
      end do
      if (mod(npaths, 2) == 1) wt(npaths) = sqrt(time)*rng_normal(rng)
    else
      do i = 1, npaths
        wt(i) = sqrt(time)*rng_normal(rng)
      end do
    end if
    xt = spot*exp((rate - 0.5_dp*sigma**2)*time + sigma*wt)

    do j = 1, size(requested)
      do i = 1, npaths
        if (present(payoff_fn)) then
          samples(i) = payoff_fn(xt(i), strike)
        else
          samples(i) = payoff_value(xt(i), strike, payoff)
        end if
      end do
      select case (trim(requested(j)))
      case ('fair_value')
      case ('delta')
        samples = samples*wt/(spot*sigma*time)
      case ('vega')
        samples = samples*(wt**2/(sigma*time) - wt - 1.0_dp/sigma)
      case ('rho')
        samples = samples*time*(wt/(sigma*time) - 1.0_dp)
      case ('theta')
        samples = -samples*(wt**2/(2.0_dp*time**2) + &
          (rate - 0.5_dp*sigma**2)*wt/(sigma*time) - &
          (1.0_dp/(2.0_dp*time) + rate))
      case ('gamma')
        samples = samples/(spot**2*sigma*time) * &
          (wt**2/(sigma*time) - wt - 1.0_dp/sigma)
      case default
        call set_error(result, greeks_invalid_argument, &
          'unsupported Malliavin European Greek')
        return
      end select
      call weighted_estimate(samples, exp(-rate*time), result, j)
    end do
  end subroutine malliavin_european_greeks

  subroutine malliavin_geometric_asian_greeks(spot, strike, rate, time, sigma, &
      dividend, payoff, requested, result, model, jump_intensity, jump_scale, &
      steps, paths, seed, antithetic, payoff_fn)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    character(len=*), intent(in) :: payoff
    character(len=*), intent(in) :: requested(:)
    type(greek_result), intent(out) :: result
    character(len=*), intent(in), optional :: model
    real(dp), intent(in), optional :: jump_intensity, jump_scale
    integer, intent(in), optional :: steps, paths, seed
    logical, intent(in), optional :: antithetic
    procedure(payoff_callback), optional :: payoff_fn
    integer :: nsteps, npaths, seed_value, i, j
    real(dp) :: lambda_value, alpha_value, dt
    logical :: anti, jump, valid
    real(dp), allocatable :: w(:, :), x(:, :), wt(:), xt(:), iw(:), i0(:)
    real(dp), allocatable :: i1(:), i2(:), ilog(:), geom(:), samples(:)

    call initialize_result(result, requested)
    nsteps = max(1, nint(time*252.0_dp)); if (present(steps)) nsteps = steps
    npaths = 10000; if (present(paths)) npaths = paths
    seed_value = 1; if (present(seed)) seed_value = seed
    anti = .false.; if (present(antithetic)) anti = antithetic
    lambda_value = 0.2_dp; if (present(jump_intensity)) lambda_value = jump_intensity
    alpha_value = 0.3_dp; if (present(jump_scale)) alpha_value = jump_scale
    jump = .false.
    if (present(model)) jump = trim(model) == 'jump_diffusion'
    call validate_mc(spot, strike, time, sigma, payoff, npaths, nsteps, result, &
      valid, present(payoff_fn))
    if (.not. valid) return
    dt = time/real(nsteps, dp)
    call simulate_paths(npaths, nsteps, time, sigma, rate - dividend, seed_value, &
      anti, jump, lambda_value, alpha_value, w, x)
    wt = w(:, nsteps + 1); xt = x(:, nsteps + 1)
    iw = calc_i(w, dt); i0 = calc_i(x, dt); i1 = calc_i_1(x, dt); i2 = calc_i_2(x, dt)
    ilog = calc_i(log(x), dt)
    allocate(geom(npaths), samples(npaths))
    geom = spot*exp(ilog/time)

    do j = 1, size(requested)
      do i = 1, npaths
        if (present(payoff_fn)) then
          samples(i) = payoff_fn(geom(i), strike)
        else
          samples(i) = payoff_value(geom(i), strike, payoff)
        end if
      end do
      select case (trim(requested(j)))
      case ('fair_value')
      case ('delta')
        samples = samples*(2.0_dp/(spot*sigma*time))*wt
      case ('rho')
        samples = samples*(wt/sigma - time)
      case ('theta')
        samples = samples*(rate + 2.0_dp*wt*ilog/(sigma*time**3) + &
          1.0_dp/time - 2.0_dp*log(xt)*wt/(sigma*time**2))
      case ('vega')
        samples = samples*(2.0_dp*wt*iw/(sigma*time**2) - 1.0_dp/sigma - wt)
      case ('gamma')
        samples = samples*(-2.0_dp*wt/(spot**2*sigma*time) + &
          4.0_dp*(wt**2 - time)/(spot**2*sigma**2*time**2))
      case default
        call set_error(result, greeks_invalid_argument, &
          'unsupported Malliavin geometric-Asian Greek')
        return
      end select
      call weighted_estimate(samples, exp(-rate*time), result, j)
    end do
  end subroutine malliavin_geometric_asian_greeks

  subroutine malliavin_asian_greeks(spot, strike, rate, time, sigma, dividend, &
      payoff, requested, result, model, jump_intensity, jump_scale, steps, &
      paths, seed, antithetic, payoff_fn, derivative_fn)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    character(len=*), intent(in) :: payoff
    character(len=*), intent(in) :: requested(:)
    type(greek_result), intent(out) :: result
    character(len=*), intent(in), optional :: model
    real(dp), intent(in), optional :: jump_intensity, jump_scale
    integer, intent(in), optional :: steps, paths, seed
    logical, intent(in), optional :: antithetic
    procedure(payoff_callback), optional :: payoff_fn, derivative_fn
    integer :: nsteps, npaths, seed_value, i, j
    real(dp) :: lambda_value, alpha_value, dt, discount
    logical :: anti, jump, valid
    real(dp), allocatable :: w(:, :), x(:, :), wt(:), xt(:), i0(:), i1(:)
    real(dp), allocatable :: i2(:), i3(:), xw(:), txw(:), average(:), samples(:)

    call initialize_result(result, requested)
    nsteps = max(1, nint(time*252.0_dp)); if (present(steps)) nsteps = steps
    npaths = 10000; if (present(paths)) npaths = paths
    seed_value = 1; if (present(seed)) seed_value = seed
    anti = .false.; if (present(antithetic)) anti = antithetic
    lambda_value = 0.2_dp; if (present(jump_intensity)) lambda_value = jump_intensity
    alpha_value = 0.3_dp; if (present(jump_scale)) alpha_value = jump_scale
    jump = .false.; if (present(model)) jump = trim(model) == 'jump_diffusion'
    call validate_mc(spot, strike, time, sigma, payoff, npaths, nsteps, result, &
      valid, present(payoff_fn))
    if (.not. valid) return
    if (present(payoff_fn) .and. .not. present(derivative_fn)) then
      do j = 1, size(requested)
        if (index(trim(requested(j)), '_d') > 0 .or. &
            trim(requested(j)) == 'gamma_kombi') then
          call set_error(result, greeks_invalid_argument, &
            'custom derivative estimators require derivative_fn')
          return
        end if
      end do
    end if

    dt = time/real(nsteps, dp); discount = exp(-rate*time)
    call simulate_paths(npaths, nsteps, time, sigma, rate - dividend, seed_value, &
      anti, jump, lambda_value, alpha_value, w, x)
    wt = w(:, nsteps + 1); xt = x(:, nsteps + 1)
    i0 = calc_i(x, dt); i1 = calc_i_1(x, dt); i2 = calc_i_2(x, dt)
    i3 = calc_i_3(x, dt); xw = calc_xw(x, w, dt); txw = calc_txw(x, w, dt)
    allocate(average(npaths), samples(npaths))
    average = spot*i0/time

    do j = 1, size(requested)
      do i = 1, npaths
        if (present(payoff_fn)) then
          samples(i) = payoff_fn(average(i), strike)
        else
          samples(i) = payoff_value(average(i), strike, payoff)
        end if
      end do
      select case (trim(requested(j)))
      case ('fair_value')
      case ('delta')
        samples = samples/(sigma*spot) * &
          (-sigma + i0*wt/i1 + sigma*i0*i2/i1**2)
      case ('delta_d')
        do i = 1, npaths
          samples(i) = derivative_value(average(i), strike, payoff, derivative_fn)*i0(i)/time
        end do
      case ('rho')
        samples = samples*(wt/sigma - time)
      case ('rho_d')
        do i = 1, npaths
          samples(i) = -time*payoff_value(average(i), strike, payoff) + &
            derivative_value(average(i), strike, payoff, derivative_fn)*spot*i1(i)/time
        end do
      case ('theta')
        samples = samples*(rate - 1.0_dp/time + &
          (i0*wt/(sigma*time) - xt*wt/sigma + time*xt)/i1 + &
          (i0*i2/time - i2*xt)/i1**2)
      case ('theta_d')
        do i = 1, npaths
          samples(i) = rate*payoff_value(average(i), strike, payoff) + &
            derivative_value(average(i), strike, payoff, derivative_fn)*spot* &
            (i0(i)/time**2 - xt(i)/time)
        end do
      case ('vega')
        samples = samples/sigma * (-(1.0_dp + sigma*wt) + &
          (wt*xw - sigma*txw)/i1 + sigma*xw*i2/i1**2)
      case ('vega_d')
        do i = 1, npaths
          samples(i) = derivative_value(average(i), strike, payoff, derivative_fn)* &
            spot*(xw(i) - sigma*i1(i))/time
        end do
      case ('gamma')
        samples = samples/(sigma**2*spot**2) * (2.0_dp*sigma**2 - &
          4.0_dp*sigma*wt*i0/i1 + &
          ((wt**2 - time)*i0 - 4.0_dp*sigma**2*i2)*i0/i1**2 + &
          sigma*(3.0_dp*wt*i2 - sigma*i3)*i0**2/i1**3 + &
          3.0_dp*sigma**2*i0**2*i2**2/i1**4)
      case ('gamma_kombi')
        do i = 1, npaths
          samples(i) = derivative_value(average(i), strike, payoff, derivative_fn) * &
            (-sigma + i0(i)*wt(i)/i1(i) + sigma*i0(i)*i2(i)/i1(i)**2) / &
            (sigma*spot)
        end do
      case default
        call set_error(result, greeks_invalid_argument, &
          'unsupported Malliavin Asian Greek')
        return
      end select
      call weighted_estimate(samples, discount, result, j)
    end do
  end subroutine malliavin_asian_greeks

  subroutine bs_malliavin_asian_greeks(spot, strike, rate, time, sigma, dividend, &
      payoff, requested, result, steps, paths, seed)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    character(len=*), intent(in) :: payoff
    character(len=*), intent(in) :: requested(:)
    type(greek_result), intent(out) :: result
    integer, intent(in), optional :: steps, paths, seed
    integer :: nsteps, npaths, seed_value, i, j, exact_index
    real(dp) :: dt, discount, exact_value
    logical :: valid
    real(dp), allocatable :: w(:, :), x(:, :), logx(:, :), wt(:), i0(:), i1(:)
    real(dp), allocatable :: i2(:), iw(:), xw(:), txw(:), geom(:), arithmetic(:)
    real(dp), allocatable :: y(:), control(:), base_payoff(:), geom_payoff(:)
    type(greek_result), allocatable :: exact
    character(len=24) :: one_name(1)

    call initialize_result(result, requested)
    nsteps = max(1, nint(time*252.0_dp)); if (present(steps)) nsteps = steps
    npaths = 1000; if (present(paths)) npaths = paths
    seed_value = 1; if (present(seed)) seed_value = seed
    call validate_mc(spot, strike, time, sigma, payoff, npaths, nsteps, result, valid)
    if (.not. valid) return
    if (trim(payoff) /= 'call' .and. trim(payoff) /= 'put') then
      call set_error(result, greeks_unknown_payoff, &
        'control-variate Asian estimator supports call or put')
      return
    end if
    dt = time/real(nsteps, dp); discount = exp(-rate*time)
    call simulate_paths(npaths, nsteps, time, sigma, rate - dividend, seed_value, &
      .false., .false., 0.0_dp, 0.0_dp, w, x)
    allocate(logx(npaths, nsteps + 1)); logx = log(x)
    wt = w(:, nsteps + 1); i0 = calc_i(x, dt); i1 = calc_i_1(x, dt)
    i2 = calc_i_2(x, dt); iw = calc_i(w, dt); xw = calc_xw(x, w, dt)
    txw = calc_txw(x, w, dt)
    geom = spot*exp(calc_i(logx, dt)/time)
    arithmetic = spot*i0/time
    allocate(y(npaths), control(npaths), base_payoff(npaths), geom_payoff(npaths))
    allocate(exact)
    do i = 1, npaths
      base_payoff(i) = payoff_value(arithmetic(i), strike, payoff)
      geom_payoff(i) = payoff_value(geom(i), strike, payoff)
    end do

    do j = 1, size(requested)
      one_name(1) = trim(requested(j))
      call bs_geometric_asian_greeks(spot, strike, rate, time, sigma, dividend, &
        payoff, one_name, exact)
      exact_index = 1
      exact_value = exact%values(exact_index)
      select case (trim(requested(j)))
      case ('fair_value')
        y = discount*base_payoff
        control = discount*geom_payoff - exact_value
      case ('delta')
        y = discount*base_payoff/(sigma*spot) * &
          (-sigma + i0*wt/i1 + sigma*i0*i2/i1**2)
        control = 2.0_dp*discount*geom_payoff*wt/(spot*sigma*time) - exact_value
      case ('rho')
        y = discount*base_payoff*(wt/sigma - time)
        control = -time*discount*geom_payoff + &
          time*discount*geom_payoff*wt/(sigma*time) - exact_value
      case ('vega')
        y = discount*base_payoff/sigma * (-(1.0_dp + sigma*wt) + &
          (wt*xw - sigma*txw)/i1 + sigma*xw*i2/i1**2)
        control = discount*geom_payoff*(2.0_dp*wt*iw/(sigma*time**2) - &
          1.0_dp/sigma - wt) - exact_value
      case default
        call set_error(result, greeks_invalid_argument, &
          'control-variate estimator supports fair_value, delta, rho, and vega')
        return
      end select
      result%values(j) = regression_intercept(y, control)
      result%standard_errors(j) = sample_stderr(y)
    end do
  end subroutine bs_malliavin_asian_greeks


  real(dp) function derivative_value(x, strike, payoff, derivative_fn) result(value)
    real(dp), intent(in) :: x, strike
    character(len=*), intent(in) :: payoff
    procedure(payoff_callback), optional :: derivative_fn
    if (present(derivative_fn)) then
      value = derivative_fn(x, strike)
    else
      value = payoff_derivative(x, strike, payoff)
    end if
  end function derivative_value

end module greeks_monte_carlo
