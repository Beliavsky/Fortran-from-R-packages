! CreditRisk calibration formulas, with sibling lbfgsb3 optimization.
! Copyright (c) 2017 Alessandro Cimarelli & Nicolò Manca. MIT license.
module creditrisk_calibration
  use r_kinds, only: dp
  use creditrisk, only: credit_curve_t, credit_ok, credit_invalid, valid_times, linear_constant, &
    black_cox, at1p, sbtv, cds
  use lbfgsb3_mod, only: lbfgsb_control_t, lbfgsb_result_t, lbfgsb_minimize_fd
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: calibration_result_t, calibrate_black_cox, calibrate_at1p, calibrate_sbtv, calibrate_cds
  public :: lbfgsb_control_t, lbfgsb_result_t

  type :: calibration_result_t
    !! Parameter order depends on model; see API.md. Error is the mean squared CDS spread error.
    real(dp), allocatable :: parameters(:)
    real(dp), allocatable :: rates(:)
    real(dp) :: error = huge(1.0_dp)
    real(dp) :: initial_error = huge(1.0_dp)
    real(dp) :: rms_error = huge(1.0_dp)
    real(dp) :: max_abs_residual = huge(1.0_dp)
    character(len=180) :: warning = 'Calibration not run: invalid inputs or starting model.' 
    integer :: status = credit_invalid
    type(lbfgsb_result_t) :: optimizer
  end type calibration_result_t

  type :: calibration_data_t
    integer :: model = 0
    real(dp) :: v0 = 1.0_dp
    real(dp) :: premium = 0.005_dp
    real(dp) :: recovery = 0.4_dp
    logical :: simplified = .false.
    real(dp) :: penalty = 1.0_dp
    real(dp), allocatable :: time(:), rate(:), quote(:), maturity(:), probability(:)
  end type calibration_data_t

contains

  function calibrate_black_cox(v0, quote, rate, time, premium, recovery, simplified, initial, control) result(fit)
    real(dp), intent(in) :: v0 !! Positive initial firm value.
    real(dp), intent(in) :: quote(:) !! Finite market CDS rates, one per maturity, in decimal units.
    real(dp), intent(in) :: rate !! Finite constant annual risk-free rate.
    real(dp), intent(in) :: time(:) !! Strictly increasing positive maturities in years.
    real(dp), intent(in), optional :: premium !! Positive CDS contractual premium; default 0.005.
    real(dp), intent(in), optional :: recovery !! CDS recovery fraction in [0,1); default 0.4.
    logical, intent(in), optional :: simplified !! True selects discrete CDS legs; false is the default.
    real(dp), intent(in), optional :: initial(:) !! Four optimizer coordinates [L,L-K,sigma,gamma]; defaults match R.
    type(lbfgsb_control_t), intent(in), optional :: control !! Optional override of all optimizer controls.
    type(calibration_result_t) :: fit
    type(calibration_data_t) :: data
    real(dp) :: start(4)
    data%model = 1
    data%v0 = v0
    data%time = time
    data%rate = [rate]
    data%quote = quote
    data%maturity = time
    call set_leg_options(data, premium, recovery, simplified)
    start = [v0 / 2.0_dp, 0.075_dp * v0, 0.2_dp, 0.08_dp]
    if (present(initial)) then
      if (size(initial) /= 4) return
      start = initial
    end if
    fit = solve(data, start, 100, control)
    ! Expose actual fitted K, correcting the R return of initial L and K.
    if (allocated(fit%parameters)) then
      if (size(fit%parameters) == 4) fit%parameters(2) = fit%parameters(1) - fit%parameters(2)
    end if
  end function calibrate_black_cox

  function calibrate_at1p(v0, quote, rate, time, premium, recovery, simplified, initial, control) result(fit)
    real(dp), intent(in) :: v0 !! Positive initial firm value, normally normalized to one.
    real(dp), intent(in) :: quote(:) !! Finite decimal market CDS spreads, one per maturity.
    real(dp), intent(in) :: rate(:) !! Finite spot rates, one per maturity.
    real(dp), intent(in) :: time(:) !! Strictly increasing positive maturities in years.
    real(dp), intent(in), optional :: premium !! Positive CDS contractual premium; default 0.005.
    real(dp), intent(in), optional :: recovery !! CDS recovery fraction in [0,1); default 0.4.
    logical, intent(in), optional :: simplified !! True selects discrete CDS legs; defaults to false.
    real(dp), intent(in), optional :: initial(:) !! Starting [H0,B,sigma(:)], length 2+size(time); R defaults [0.4,0.02,0.2,...].
    type(lbfgsb_control_t), intent(in), optional :: control !! Optional override of all optimizer controls.
    type(calibration_result_t) :: fit
    type(calibration_data_t) :: data
    real(dp), allocatable :: start(:)
    data%model = 2
    data%v0 = v0
    data%time = time
    data%rate = rate
    data%quote = quote
    data%maturity = time
    call set_leg_options(data, premium, recovery, simplified)
    allocate(start(2 + size(time)))
    start = 0.2_dp
    start(1:2) = [0.4_dp, 0.02_dp]
    if (present(initial)) then
      if (size(initial) /= size(start)) return
      start = initial
    end if
    fit = solve(data, start, 5000, control)
  end function calibrate_at1p

  function calibrate_sbtv(v0, probability, quote, rate, time, initial, control) result(fit)
    real(dp), intent(in) :: v0 !! Positive initial firm value.
    real(dp), intent(in) :: probability(:) !! Nonnegative fixed scenario weights summing to one within 1e-12.
    real(dp), intent(in) :: quote(:) !! Finite decimal CDS market spreads, one per maturity.
    real(dp), intent(in) :: rate(:) !! Finite spot rates, one per maturity.
    real(dp), intent(in) :: time(:) !! Strictly increasing positive maturities in years.
    real(dp), intent(in), optional :: initial(:) !! Starting [B,sigma(:),H(:)], length 1+dates+scenarios; R defaults.
    type(lbfgsb_control_t), intent(in), optional :: control !! Optional override of all optimizer controls.
    type(calibration_result_t) :: fit
    type(calibration_data_t) :: data
    real(dp), allocatable :: start(:)
    integer :: n, k, i
    data%model = 3
    data%v0 = v0
    data%time = time
    data%rate = rate
    data%quote = quote
    data%maturity = time
    data%probability = probability
    n = size(time)
    k = size(probability)
    if (k == 0) return
    allocate(start(1 + n + k))
    start(1) = 0.02_dp
    start(2:1 + n) = 0.2_dp
    do i = 1, k
      start(1 + n + i) = 0.4_dp * v0 * 1.5_dp ** (i - 1)
    end do
    if (present(initial)) then
      if (size(initial) /= size(start)) return
      start = initial
    end if
    ! R ignores ... in this calibrator; retain default CDS leg settings.
    fit = solve(data, start, 5000, control)
  end function calibrate_sbtv

  function calibrate_cds(rate, time, maturity, quote, premium, recovery, simplified, initial, control) result(fit)
    real(dp), intent(in) :: rate(:) !! Finite spot rates on payment dates; scalar array broadcasts.
    real(dp), intent(in) :: time(:) !! Strictly increasing positive premium-payment dates in years.
    real(dp), intent(in) :: maturity(:) !! Strictly increasing positive market-quote dates and hazard knots.
    real(dp), intent(in) :: quote(:) !! Finite decimal market spreads, one per quote maturity.
    real(dp), intent(in), optional :: premium !! Positive CDS contractual premium; default 0.005.
    real(dp), intent(in), optional :: recovery !! CDS recovery fraction in [0,1); default 0.4.
    logical, intent(in), optional :: simplified !! True selects discrete CDS legs; defaults to false.
    real(dp), intent(in), optional :: initial(:) !! Starting hazards on maturity knots, default 0.1 at every knot.
    type(lbfgsb_control_t), intent(in), optional :: control !! Optional override of all optimizer controls.
    type(calibration_result_t) :: fit
    type(calibration_data_t) :: data
    real(dp), allocatable :: start(:)
    data%model = 4
    data%time = time
    data%rate = rate
    data%quote = quote
    data%maturity = maturity
    call set_leg_options(data, premium, recovery, simplified)
    allocate(start(size(maturity)))
    start = 0.1_dp
    if (present(initial)) then
      if (size(initial) /= size(start)) return
      start = initial
    end if
    fit = solve(data, start, 5000, control)
  end function calibrate_cds

  pure subroutine set_leg_options(data, premium, recovery, simplified)
    type(calibration_data_t), intent(inout) :: data !! Calibration context receiving optional CDS-leg settings.
    real(dp), intent(in), optional :: premium !! Positive contractual premium replacing the default when present.
    real(dp), intent(in), optional :: recovery !! Recovery fraction in [0,1), replacing the default when present.
    logical, intent(in), optional :: simplified !! Discrete-leg switch replacing the default when present.
    if (present(premium)) data%premium = premium
    if (present(recovery)) data%recovery = recovery
    if (present(simplified)) data%simplified = simplified
  end subroutine set_leg_options

  function solve(data, start, max_iterations, control) result(fit)
    type(calibration_data_t), intent(inout) :: data !! Validated market data and model-specific objective settings.
    real(dp), intent(in) :: start(:) !! Nonnegative starting optimizer coordinates.
    integer, intent(in) :: max_iterations !! R-compatible default iteration cap, overridden by control.
    type(lbfgsb_control_t), intent(in), optional :: control !! Optional complete replacement of optimizer settings.
    type(calibration_result_t) :: fit
    type(lbfgsb_control_t) :: settings
    logical :: valid
    if (.not. valid_times(data%time) .or. .not. valid_times(data%maturity)) return
    if (size(data%quote) /= size(data%maturity) .or. size(start) == 0) return
    if (.not. all(ieee_is_finite(data%quote)) .or. .not. all(ieee_is_finite(data%rate))) return
    if (.not. all(ieee_is_finite(start)) .or. any(start < 0.0_dp)) return
    if (data%model == 2 .or. data%model == 3) then
      if (size(data%rate) /= size(data%time)) return
    else
      if (size(data%rate) /= 1 .and. size(data%rate) /= size(data%time)) return
    end if
    fit%parameters = start
    call model_rates(start, data, fit%rates, valid)
    if (.not. valid) return
    fit%initial_error = mean_squared_error(fit%rates, data%quote)
    if (.not. ieee_is_finite(fit%initial_error)) return
    ! A huge artificial penalty can collapse line-search steps to zero and
    ! falsely satisfy relative-function convergence. Scale to the starting MSE.
    data%penalty = max(1.0_dp, 100.0_dp * fit%initial_error)
    settings%memory = 5
    settings%factr = 1.0e7_dp
    settings%pgtol = 0.0_dp
    settings%abstol = 0.0_dp
    settings%reltol = 0.0_dp
    settings%finite_difference_step = 0.001_dp
    settings%max_iterations = max_iterations
    settings%max_evaluations = huge(0)
    if (present(control)) settings = control
    call lbfgsb_minimize_fd(fit%parameters, objective, fit%optimizer, lower=[0.0_dp], control=settings, user_data=data)
    call model_rates(fit%parameters, data, fit%rates, valid)
    if (.not. valid) return
    fit%error = mean_squared_error(fit%rates, data%quote)
    if (.not. ieee_is_finite(fit%error)) return
    fit%rms_error = sqrt(fit%error)
    fit%max_abs_residual = maxval(abs(fit%rates - data%quote))
    fit%warning = ''
    if (.not. fit%optimizer%success) fit%warning = 'Optimizer did not converge. Inspect optimizer diagnostics.'
    if (fit%max_abs_residual > 0.0001_dp) then
      fit%warning = trim(fit%warning) // ' Maximum quote residual exceeds one basis point. Inspect fit quality.'
      fit%warning = adjustl(fit%warning)
    end if
    fit%status = credit_ok
  end function solve

  pure real(dp) function mean_squared_error(rates, quotes) result(error)
    real(dp), intent(in) :: rates(:) !! Finite fitted spreads with the same length as quotes.
    real(dp), intent(in) :: quotes(:) !! Finite market spreads defining a nonempty objective.
    real(dp) :: scale
    error = ieee_value(0.0_dp, ieee_quiet_nan)
    if (size(rates) == 0 .or. size(rates) /= size(quotes)) return
    if (.not. all(ieee_is_finite(rates)) .or. .not. all(ieee_is_finite(quotes))) return
    ! Bound before subtraction/squaring, leaving room for the 100-fold penalty.
    if (max(maxval(abs(rates)), maxval(abs(quotes))) > sqrt(huge(1.0_dp) / 400.0_dp)) return
    scale = maxval(abs(rates - quotes))
    error = 0.0_dp
    if (scale == 0.0_dp) return
    error = scale ** 2 * (sum(((rates - quotes) / scale) ** 2) / real(size(quotes), dp))
  end function mean_squared_error

  function objective(parameters, user_data) result(error)
    real(dp), intent(in) :: parameters(:) !! Nonnegative L-BFGS-B trial coordinates for the selected model.
    class(*), intent(inout), optional :: user_data !! Required calibration_data_t market context passed by the optimizer.
    real(dp) :: error
    real(dp), allocatable :: rates(:)
    logical :: valid
    error = 1.0_dp
    if (.not. present(user_data)) return
    select type (user_data)
    type is (calibration_data_t)
      error = user_data%penalty
      call model_rates(parameters, user_data, rates, valid)
      if (valid) error = mean_squared_error(rates, user_data%quote)
      if (.not. ieee_is_finite(error)) error = user_data%penalty
    end select
  end function objective

  pure subroutine model_rates(parameters, data, rates, valid)
    real(dp), intent(in) :: parameters(:) !! Trial coordinates in the documented model-specific optimizer order.
    type(calibration_data_t), intent(in) :: data !! Market dates, rates, fixed parameters and CDS-leg settings.
    real(dp), allocatable, intent(out) :: rates(:) !! Fitted CDS rates on quote dates; empty when invalid.
    logical, intent(out) :: valid !! True only for a valid model and finite fitted rates.
    type(credit_curve_t) :: structural, legs
    real(dp), allocatable :: intensity(:)
    integer :: n
    valid = .false.
    allocate(rates(0))
    n = size(data%time)
    select case (data%model)
    case (1)
      structural = black_cox(parameters(1), data%v0, parameters(3), data%rate(1), parameters(4), &
        data%time, parameters(1) - parameters(2))
      if (structural%status /= credit_ok) return
      intensity = structural%values(:,5)
    case (2)
      structural = at1p(data%v0, parameters(1), parameters(2), parameters(3:), data%rate, data%time)
      if (structural%status /= credit_ok) return
      intensity = structural%values(:,5)
    case (3)
      structural = sbtv(data%v0, parameters(2 + n:), data%probability, parameters(1), &
        parameters(2:1 + n), data%rate, data%time)
      if (structural%status /= credit_ok) return
      intensity = structural%values(:,3)
    case (4)
      intensity = linear_constant(data%maturity, parameters, data%time)
    case default
      return
    end select
    legs = cds(data%time, intensity, data%rate, data%premium, data%recovery, data%simplified)
    if (legs%status /= credit_ok) return
    if (data%model == 4) then
      ! calibrate.cds interpolates rates, whereas cds2 takes a ratio of interpolated legs.
      rates = linear_constant(data%time, legs%values(:,5), data%maturity)
    else
      rates = legs%values(:,5)
    end if
    valid = all(ieee_is_finite(rates))
  end subroutine model_rates

end module creditrisk_calibration
