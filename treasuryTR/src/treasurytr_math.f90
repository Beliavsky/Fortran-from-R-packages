! SPDX-License-Identifier: MIT
! Copyright (c) 2021 treasuryTR authors

module treasurytr_math
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
  use treasurytr_kinds, only : dp
  implicit none
  private

  public :: convexity
  public :: mod_duration
  public :: period_total_return

contains

  pure elemental function mod_duration(yield_rate, maturity, source_compatible) result(value)
    real(dp), intent(in) :: yield_rate
    real(dp), intent(in) :: maturity
    logical, intent(in), optional :: source_compatible
    real(dp) :: value
    logical :: use_source
    real(dp) :: n

    use_source = .true.
    if (present(source_compatible)) use_source = source_compatible

    if (.not. ieee_is_finite(yield_rate) .or. .not. ieee_is_finite(maturity)) then
      value = quiet_nan()
      return
    end if
    if (maturity <= 0.0_dp .or. yield_rate <= -2.0_dp) then
      value = quiet_nan()
      return
    end if

    if (yield_rate == 0.0_dp) then
      if (use_source) then
        value = quiet_nan()
      else
        value = maturity
      end if
      return
    end if

    if (use_source) then
      value = (1.0_dp - 1.0_dp / (1.0_dp + 0.5_dp * yield_rate) ** &
        (2.0_dp * maturity)) / yield_rate
    else if (abs(yield_rate) < 1.0e-5_dp) then
      n = 2.0_dp * maturity
      value = 0.5_dp * n - n * (n + 1.0_dp) * yield_rate / 8.0_dp + &
        n * (n + 1.0_dp) * (n + 2.0_dp) * yield_rate ** 2 / 48.0_dp - &
        n * (n + 1.0_dp) * (n + 2.0_dp) * (n + 3.0_dp) * &
        yield_rate ** 3 / 384.0_dp
    else
      value = -expm1_safe(-2.0_dp * maturity * log1p_safe(0.5_dp * yield_rate)) / &
        yield_rate
    end if
  end function mod_duration

  pure elemental function convexity(yield_rate, maturity, source_compatible) result(value)
    real(dp), intent(in) :: yield_rate
    real(dp), intent(in) :: maturity
    logical, intent(in), optional :: source_compatible
    real(dp) :: value
    logical :: use_source
    real(dp) :: n, z, log_z

    use_source = .true.
    if (present(source_compatible)) use_source = source_compatible

    if (.not. ieee_is_finite(yield_rate) .or. .not. ieee_is_finite(maturity)) then
      value = quiet_nan()
      return
    end if
    if (maturity <= 0.0_dp .or. yield_rate <= -2.0_dp) then
      value = quiet_nan()
      return
    end if

    n = 2.0_dp * maturity
    if (yield_rate == 0.0_dp) then
      if (use_source) then
        value = quiet_nan()
      else
        value = n * (n + 1.0_dp) / 4.0_dp
      end if
      return
    end if

    if (use_source) then
      z = 1.0_dp + 0.5_dp * yield_rate
      value = 2.0_dp / yield_rate ** 2 * (1.0_dp - z ** (-n)) - &
        n / yield_rate * z ** (-n - 1.0_dp)
    else if (abs(yield_rate) < 1.0e-4_dp) then
      value = n * (n + 1.0_dp) / 4.0_dp - &
        n * (n + 1.0_dp) * (n + 2.0_dp) * yield_rate / 12.0_dp + &
        n * (n + 1.0_dp) * (n + 2.0_dp) * (n + 3.0_dp) * &
        yield_rate ** 2 / 64.0_dp - &
        n * (n + 1.0_dp) * (n + 2.0_dp) * (n + 3.0_dp) * &
        (n + 4.0_dp) * yield_rate ** 3 / 480.0_dp
    else
      log_z = log1p_safe(0.5_dp * yield_rate)
      value = 2.0_dp / yield_rate ** 2 * &
        (-expm1_safe(-n * log_z)) - n / yield_rate * exp(-(n + 1.0_dp) * log_z)
    end if
  end function convexity

  pure elemental function period_total_return(current_yield, previous_yield, maturity, &
      scale, mdur_current, convex_current, source_compatible) result(value)
    real(dp), intent(in) :: current_yield
    real(dp), intent(in) :: previous_yield
    real(dp), intent(in) :: maturity
    real(dp), intent(in) :: scale
    real(dp), intent(in), optional :: mdur_current
    real(dp), intent(in), optional :: convex_current
    logical, intent(in), optional :: source_compatible
    real(dp) :: value
    real(dp) :: delta_y, duration_value, convexity_value, income
    logical :: use_source

    use_source = .true.
    if (present(source_compatible)) use_source = source_compatible

    if (.not. ieee_is_finite(current_yield) .or. &
        .not. ieee_is_finite(previous_yield) .or. &
        .not. ieee_is_finite(maturity) .or. .not. ieee_is_finite(scale)) then
      value = quiet_nan()
      return
    end if
    if (maturity <= 0.0_dp .or. scale <= 0.0_dp .or. previous_yield <= -1.0_dp) then
      value = quiet_nan()
      return
    end if

    if (present(mdur_current)) then
      duration_value = mdur_current
    else
      duration_value = mod_duration(current_yield, maturity, use_source)
    end if
    if (present(convex_current)) then
      convexity_value = convex_current
    else
      convexity_value = convexity(current_yield, maturity, use_source)
    end if

    if (.not. ieee_is_finite(duration_value) .or. &
        .not. ieee_is_finite(convexity_value)) then
      value = quiet_nan()
      return
    end if

    delta_y = current_yield - previous_yield
    if (use_source) then
      income = (1.0_dp + previous_yield) ** (1.0_dp / scale) - 1.0_dp
    else
      income = expm1_safe(log1p_safe(previous_yield) / scale)
    end if
    value = -duration_value * delta_y + 0.5_dp * convexity_value * delta_y ** 2 + income
  end function period_total_return

  pure elemental function log1p_safe(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value

    if (abs(x) < 1.0e-8_dp) then
      value = x - 0.5_dp * x ** 2 + x ** 3 / 3.0_dp - 0.25_dp * x ** 4 + &
        0.2_dp * x ** 5
    else
      value = log(1.0_dp + x)
    end if
  end function log1p_safe

  pure elemental function expm1_safe(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value

    if (abs(x) < 1.0e-6_dp) then
      value = x + 0.5_dp * x ** 2 + x ** 3 / 6.0_dp + x ** 4 / 24.0_dp + &
        x ** 5 / 120.0_dp
    else
      value = exp(x) - 1.0_dp
    end if
  end function expm1_safe

  pure function quiet_nan() result(value)
    real(dp) :: value

    value = ieee_value(0.0_dp, ieee_quiet_nan)
  end function quiet_nan

end module treasurytr_math
