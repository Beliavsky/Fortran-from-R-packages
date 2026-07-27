! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
module greeks_malliavin
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use greeks_kinds, only: dp
  use greeks_types, only: mc_greeks_result, greeks_result, payoff_call, payoff_put, &
    payoff_cash_call, payoff_cash_put, payoff_asset_call, payoff_asset_put, &
    model_black_scholes, model_jump_diffusion, payoff_callback, &
    payoff_derivative_callback, jump_sampler_callback
  use greeks_rng, only: seed_rng, normal_random
  use greeks_math, only: sample_mean, sample_se, control_variate_stats
  use greeks_paths, only: simulate_paths, trapezoid_integral, weighted_integral, &
    integral_xw, integral_txw
  use greeks_geometric_asian, only: bs_geometric_asian_greeks
  implicit none
  private
  public :: malliavin_european_greeks
  public :: malliavin_geometric_asian_greeks
  public :: malliavin_asian_greeks
  public :: bs_malliavin_asian_greeks
contains
  pure function builtin_payoff(x, strike, payoff) result(value)
    real(dp), intent(in) :: x, strike
    integer, intent(in) :: payoff
    real(dp) :: value
    select case (payoff)
    case (payoff_call)
      value = max(x-strike,0.0_dp)
    case (payoff_put)
      value = max(strike-x,0.0_dp)
    case (payoff_cash_call)
      if (x >= strike) then
        value = 1.0_dp
      else
        value = 0.0_dp
      end if
    case (payoff_cash_put)
      if (x <= strike) then
        value = 1.0_dp
      else
        value = 0.0_dp
      end if
    case (payoff_asset_call)
      if (x >= strike) then
        value = x
      else
        value = 0.0_dp
      end if
    case (payoff_asset_put)
      if (x <= strike) then
        value = x
      else
        value = 0.0_dp
      end if
    case default
      value = 0.0_dp
    end select
  end function builtin_payoff

  pure function builtin_dpayoff(x, strike, payoff) result(value)
    real(dp), intent(in) :: x, strike
    integer, intent(in) :: payoff
    real(dp) :: value
    select case (payoff)
    case (payoff_call)
      if (x > strike) then
        value = 1.0_dp
      else
        value = 0.0_dp
      end if
    case (payoff_put)
      if (x < strike) then
        value = -1.0_dp
      else
        value = 0.0_dp
      end if
    case default
      value = 0.0_dp
    end select
  end function builtin_dpayoff

  subroutine payoff_vector(x, strike, payoff, values, payoff_fn)
    real(dp), intent(in) :: x(:), strike
    integer, intent(in) :: payoff
    real(dp), intent(out) :: values(size(x))
    procedure(payoff_callback), optional :: payoff_fn
    integer :: i
    if (present(payoff_fn)) then
      do i = 1, size(x)
        values(i) = payoff_fn(x(i),strike)
      end do
    else
      do i = 1, size(x)
        values(i) = builtin_payoff(x(i),strike,payoff)
      end do
    end if
  end subroutine payoff_vector

  subroutine dpayoff_vector(x, strike, payoff, values, dpayoff_fn)
    real(dp), intent(in) :: x(:), strike
    integer, intent(in) :: payoff
    real(dp), intent(out) :: values(size(x))
    procedure(payoff_derivative_callback), optional :: dpayoff_fn
    integer :: i
    if (present(dpayoff_fn)) then
      do i = 1, size(x)
        values(i) = dpayoff_fn(x(i),strike)
      end do
    else
      do i = 1, size(x)
        values(i) = builtin_dpayoff(x(i),strike,payoff)
      end do
    end if
  end subroutine dpayoff_vector

  subroutine set_mc_component(contribution, estimate, standard_error)
    real(dp), intent(in) :: contribution(:)
    real(dp), intent(out) :: estimate, standard_error
    estimate = sample_mean(contribution)
    standard_error = sample_se(contribution)
  end subroutine set_mc_component

  function malliavin_european_greeks(spot, strike, rate, time, sigma, payoff, paths, seed, &
      antithetic, payoff_fn) result(res)
    real(dp), intent(in) :: spot, strike, rate, time, sigma
    integer, intent(in) :: payoff
    integer, intent(in), optional :: paths, seed
    logical, intent(in), optional :: antithetic
    procedure(payoff_callback), optional :: payoff_fn
    type(mc_greeks_result) :: res
    integer :: n, i, half, source
    logical :: anti
    real(dp), allocatable :: w(:), terminal(:), pval(:), contribution(:)
    real(dp) :: disc

    n = 10000
    if (present(paths)) n = paths
    anti = .false.
    if (present(antithetic)) anti = antithetic
    if (n < 2 .or. spot <= 0.0_dp .or. strike <= 0.0_dp .or. time <= 0.0_dp .or. sigma <= 0.0_dp) then
      res%ok = .false.
      res%message = 'invalid dimensions or Black-Scholes parameters'
      res%estimate%ok = .false.
      res%estimate%message = res%message
      return
    end if
    if (present(seed)) then
      call seed_rng(seed)
    else
      call seed_rng(1)
    end if
    allocate(w(n),terminal(n),pval(n),contribution(n))
    half = n
    if (anti) half = (n+1)/2
    do i = 1, half
      w(i) = sqrt(time)*normal_random()
    end do
    if (anti) then
      do i = half+1, n
        source = i-half
        w(i) = -w(source)
      end do
    end if
    terminal = spot*exp((rate-0.5_dp*sigma*sigma)*time+sigma*w)
    call payoff_vector(terminal,strike,payoff,pval,payoff_fn)
    disc = exp(-rate*time)
    contribution = disc*pval
    call set_mc_component(contribution,res%estimate%fair_value,res%standard_error%fair_value)
    contribution = disc*pval*w/(spot*sigma*time)
    call set_mc_component(contribution,res%estimate%delta,res%standard_error%delta)
    contribution = disc*pval*(w*w/(sigma*time)-w-1.0_dp/sigma)
    call set_mc_component(contribution,res%estimate%vega,res%standard_error%vega)
    contribution = disc*pval*(w/sigma-time)
    call set_mc_component(contribution,res%estimate%rho,res%standard_error%rho)
    contribution = -disc*pval*(w*w/(2.0_dp*time*time) + &
      (rate-0.5_dp*sigma*sigma)*w/(sigma*time) - (1.0_dp/(2.0_dp*time)+rate))
    call set_mc_component(contribution,res%estimate%theta,res%standard_error%theta)
    contribution = disc*pval*(w*w/(sigma*time)-w-1.0_dp/sigma)/(spot*spot*sigma*time)
    call set_mc_component(contribution,res%estimate%gamma,res%standard_error%gamma)
    res%paths = n
  end function malliavin_european_greeks

  function malliavin_geometric_asian_greeks(spot, strike, rate, time, sigma, dividend, payoff, &
      steps, paths, seed, antithetic, model, jump_lambda, jump_scale, payoff_fn, jump_sampler) result(res)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    integer, intent(in) :: payoff
    integer, intent(in), optional :: steps, paths, seed, model
    logical, intent(in), optional :: antithetic
    real(dp), intent(in), optional :: jump_lambda, jump_scale
    procedure(payoff_callback), optional :: payoff_fn
    procedure(jump_sampler_callback), optional :: jump_sampler
    type(mc_greeks_result) :: res
    integer :: ns, np, model_id
    logical :: anti
    real(dp) :: lambda, alpha, dt, disc
    real(dp), allocatable :: w(:,:), x(:,:), log_sx(:,:), geom(:), pval(:)
    real(dp), allocatable :: iw(:), ilog(:), contribution(:)

    ns = max(1,nint(time*252.0_dp))
    if (present(steps)) ns = steps
    np = 10000
    if (present(paths)) np = paths
    anti = .false.
    if (present(antithetic)) anti = antithetic
    model_id = model_black_scholes
    if (present(model)) model_id = model
    lambda = 0.2_dp
    if (present(jump_lambda)) lambda = jump_lambda
    alpha = 0.3_dp
    if (present(jump_scale)) alpha = jump_scale
    if (ns < 1 .or. np < 2 .or. spot <= 0.0_dp .or. strike <= 0.0_dp .or. &
        time <= 0.0_dp .or. sigma <= 0.0_dp) then
      res%ok = .false.
      res%message = 'invalid dimensions or model parameters'
      return
    end if
    if (present(seed)) then
      call seed_rng(seed)
    else
      call seed_rng(1)
    end if
    if (present(jump_sampler)) then
      call simulate_paths(np,ns,time,sigma,rate-dividend,model_id,lambda,alpha,anti,w,x,jump_sampler)
    else
      call simulate_paths(np,ns,time,sigma,rate-dividend,model_id,lambda,alpha,anti,w,x)
    end if
    dt = time/real(ns,dp)
    allocate(log_sx(np,0:ns),geom(np),pval(np),contribution(np))
    log_sx = log(spot*x)
    ilog = trapezoid_integral(log_sx,dt)
    geom = exp(ilog/time)
    call payoff_vector(geom,strike,payoff,pval,payoff_fn)
    iw = trapezoid_integral(w,dt)
    disc = exp(-rate*time)
    contribution = disc*pval
    call set_mc_component(contribution,res%estimate%fair_value,res%standard_error%fair_value)
    contribution = disc*pval*(2.0_dp/(spot*sigma*time))*w(:,ns)
    call set_mc_component(contribution,res%estimate%delta,res%standard_error%delta)
    contribution = disc*pval*(w(:,ns)/sigma-time)
    call set_mc_component(contribution,res%estimate%rho,res%standard_error%rho)
    contribution = disc*pval*(2.0_dp*w(:,ns)*iw/(sigma*time*time)-1.0_dp/sigma-w(:,ns))
    call set_mc_component(contribution,res%estimate%vega,res%standard_error%vega)
    contribution = disc*pval*(-2.0_dp*w(:,ns)/(spot*spot*sigma*time) + &
      4.0_dp*(w(:,ns)*w(:,ns)-time)/(spot*spot*sigma*sigma*time))
    call set_mc_component(contribution,res%estimate%gamma,res%standard_error%gamma)
    contribution = disc*pval*(rate + 2.0_dp*w(:,ns)*trapezoid_integral(log(x),dt)/(sigma*time**3) - &
      1.0_dp/time - 2.0_dp*log(x(:,ns))*w(:,ns)/(sigma*time*time) + 2.0_dp/time)
    call set_mc_component(contribution,res%estimate%theta,res%standard_error%theta)
    res%paths = np
    res%steps = ns
  end function malliavin_geometric_asian_greeks

  function malliavin_asian_greeks(spot, strike, rate, time, sigma, dividend, payoff, &
      steps, paths, seed, antithetic, model, jump_lambda, jump_scale, payoff_fn, dpayoff_fn, &
      jump_sampler) result(res)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    integer, intent(in) :: payoff
    integer, intent(in), optional :: steps, paths, seed, model
    logical, intent(in), optional :: antithetic
    real(dp), intent(in), optional :: jump_lambda, jump_scale
    procedure(payoff_callback), optional :: payoff_fn
    procedure(payoff_derivative_callback), optional :: dpayoff_fn
    procedure(jump_sampler_callback), optional :: jump_sampler
    type(mc_greeks_result) :: res
    integer :: ns, np, model_id
    logical :: anti
    real(dp) :: lambda, alpha, dt, disc
    real(dp), allocatable :: w(:,:), x(:,:), avg(:), pval(:), dpval(:), c(:)
    real(dp), allocatable :: i0(:), i1(:), i2(:), i3(:), xw(:), txw(:)

    ns = max(1,nint(time*252.0_dp)); if (present(steps)) ns = steps
    np = 10000; if (present(paths)) np = paths
    anti = .false.; if (present(antithetic)) anti = antithetic
    model_id = model_black_scholes; if (present(model)) model_id = model
    lambda = 0.2_dp; if (present(jump_lambda)) lambda = jump_lambda
    alpha = 0.3_dp; if (present(jump_scale)) alpha = jump_scale
    if (ns < 1 .or. np < 2 .or. spot <= 0.0_dp .or. strike <= 0.0_dp .or. &
        time <= 0.0_dp .or. sigma <= 0.0_dp) then
      res%ok = .false.; res%message = 'invalid dimensions or model parameters'; return
    end if
    if (present(seed)) then; call seed_rng(seed); else; call seed_rng(1); end if
    if (present(jump_sampler)) then
      call simulate_paths(np,ns,time,sigma,rate-dividend,model_id,lambda,alpha,anti,w,x,jump_sampler)
    else
      call simulate_paths(np,ns,time,sigma,rate-dividend,model_id,lambda,alpha,anti,w,x)
    end if
    dt = time/real(ns,dp)
    i0 = trapezoid_integral(x,dt)
    i1 = weighted_integral(x,dt,1)
    i2 = weighted_integral(x,dt,2)
    i3 = weighted_integral(x,dt,3)
    xw = integral_xw(x,w,dt)
    txw = integral_txw(x,w,dt)
    allocate(avg(np),pval(np),dpval(np),c(np))
    avg = spot*i0/time
    call payoff_vector(avg,strike,payoff,pval,payoff_fn)
    call dpayoff_vector(avg,strike,payoff,dpval,dpayoff_fn)
    disc = exp(-rate*time)
    c = disc*pval
    call set_mc_component(c,res%estimate%fair_value,res%standard_error%fair_value)
    c = disc*pval/(sigma*spot)*(-sigma+i0*w(:,ns)/i1+sigma*i0*i2/(i1*i1))
    call set_mc_component(c,res%estimate%delta,res%standard_error%delta)
    c = disc*dpval*i0/time
    call set_mc_component(c,res%estimate%delta_d,res%standard_error%delta_d)
    c = disc*pval*(w(:,ns)/sigma-time)
    call set_mc_component(c,res%estimate%rho,res%standard_error%rho)
    c = -time*disc*pval + disc*dpval*spot*i1/time
    call set_mc_component(c,res%estimate%rho_d,res%standard_error%rho_d)
    c = disc*pval*(rate-1.0_dp/time + ((i0*w(:,ns)/(sigma*time)-x(:,ns)*w(:,ns)/sigma + &
      time*x(:,ns))/i1) + (i0*i2/time-i2*x(:,ns))/(i1*i1))
    call set_mc_component(c,res%estimate%theta,res%standard_error%theta)
    c = rate*disc*pval + disc*dpval*spot*(i0/time**2-x(:,ns)/time)
    call set_mc_component(c,res%estimate%theta_d,res%standard_error%theta_d)
    c = disc*pval/sigma*(-(1.0_dp+sigma*w(:,ns)) + &
      (w(:,ns)*xw-sigma*txw)/i1 + sigma*xw*i2/(i1*i1))
    call set_mc_component(c,res%estimate%vega,res%standard_error%vega)
    c = disc*dpval*(spot/time)*(xw-sigma*i1)
    call set_mc_component(c,res%estimate%vega_d,res%standard_error%vega_d)
    c = disc*pval/(sigma*sigma*spot*spot)*(2.0_dp*sigma*sigma - &
      4.0_dp*sigma*w(:,ns)*i0/i1 + (((w(:,ns)**2-time)*i0-4.0_dp*sigma*sigma*i2)*i0)/(i1*i1) + &
      sigma*(3.0_dp*w(:,ns)*i2-sigma*i3)*i0*i0/(i1**3) + &
      3.0_dp*sigma*sigma*i0*i0*i2*i2/(i1**4))
    call set_mc_component(c,res%estimate%gamma,res%standard_error%gamma)
    c = disc*dpval/(sigma*spot)*(-sigma+i0*w(:,ns)/i1+sigma*i0*i2/(i1*i1))
    call set_mc_component(c,res%estimate%gamma_kombi,res%standard_error%gamma_kombi)
    res%paths=np; res%steps=ns
  end function malliavin_asian_greeks

  function bs_malliavin_asian_greeks(spot, strike, rate, time, sigma, dividend, payoff, &
      steps, paths, seed) result(res)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    integer, intent(in) :: payoff
    integer, intent(in), optional :: steps, paths, seed
    type(mc_greeks_result) :: res
    integer :: ns, np
    real(dp) :: dt, disc
    real(dp), allocatable :: w(:,:), x(:,:), logx(:,:), i0(:), i1(:), i2(:)
    real(dp), allocatable :: iw(:), xw(:), txw(:), avg(:), geom(:), p_ar(:), p_geom(:)
    real(dp), allocatable :: y(:), cv(:)
    type(greeks_result) :: exact

    ns=max(1,nint(time*252.0_dp)); if (present(steps)) ns=steps
    np=1000; if (present(paths)) np=paths
    if (present(seed)) then; call seed_rng(seed); else; call seed_rng(1); end if
    call simulate_paths(np,ns,time,sigma,rate-dividend,model_black_scholes,0.0_dp,0.0_dp,.false.,w,x)
    dt=time/real(ns,dp); disc=exp(-rate*time)
    i0=trapezoid_integral(x,dt); i1=weighted_integral(x,dt,1); i2=weighted_integral(x,dt,2)
    iw=trapezoid_integral(w,dt); xw=integral_xw(x,w,dt); txw=integral_txw(x,w,dt)
    allocate(logx(np,0:ns),avg(np),geom(np),p_ar(np),p_geom(np),y(np),cv(np))
    logx=log(x); avg=spot*i0/time; geom=spot*exp(trapezoid_integral(logx,dt)/time)
    call payoff_vector(avg,strike,payoff,p_ar)
    call payoff_vector(geom,strike,payoff,p_geom)
    exact=bs_geometric_asian_greeks(spot,strike,rate,time,sigma,dividend,payoff)
    y=disc*p_ar; cv=disc*p_geom-exact%fair_value
    call control_variate_stats(y,cv,res%estimate%fair_value,res%standard_error%fair_value)
    y=disc*p_ar/(sigma*spot)*(-sigma+i0*w(:,ns)/i1+sigma*i0*i2/(i1*i1))
    cv=2.0_dp*disc*p_geom*w(:,ns)/(spot*sigma*time)-exact%delta
    call control_variate_stats(y,cv,res%estimate%delta,res%standard_error%delta)
    y=disc*p_ar*(w(:,ns)/sigma-time)
    cv=(-time*disc*p_geom+spot*time/2.0_dp*(2.0_dp*disc*p_geom*w(:,ns)/(spot*sigma*time)))-exact%rho
    call control_variate_stats(y,cv,res%estimate%rho,res%standard_error%rho)
    y=disc*p_ar/sigma*(-(1.0_dp+sigma*w(:,ns)) + &
      (w(:,ns)*xw-sigma*txw)/i1+sigma*xw*i2/(i1*i1))
    cv=disc*p_geom*(2.0_dp*w(:,ns)*iw/(sigma*time*time)-1.0_dp/sigma-w(:,ns))-exact%vega
    call control_variate_stats(y,cv,res%estimate%vega,res%standard_error%vega)
    res%paths=np; res%steps=ns
  end function bs_malliavin_asian_greeks
end module greeks_malliavin
