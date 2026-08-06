! SPDX-License-Identifier: GPL-3.0-or-later
module robstattm_psi
  use robstattm_kinds, only : dp, pi
  use robstattm_utils, only : lower_string, median_absolute, mean_value
  implicit none
  private
  public :: rho_value, rho_prime, rho_second, rho_weight
  public :: tuning_huber, tuning_bisquare, tuning_opt, tuning_mopt
  public :: tuning_for_efficiency, scale_m, inverse_robust_r_squared
  public :: gaussian_efficiency, rho_shr, weight_shr
contains
  elemental function normal_density(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    value = exp(-0.5_dp * x * x) / sqrt(2.0_dp * pi)
  end function normal_density

  elemental function rho_value(u, family, cc, standardize) result(value)
    real(dp), intent(in) :: u, cc
    character(len=*), intent(in) :: family
    logical, intent(in), optional :: standardize
    real(dp) :: value, z, raw
    logical :: standardized
    character(len=len(family)) :: name
    name = lower_string(trim(adjustl(family)))
    standardized = .true.
    if (present(standardize)) standardized = standardize
    select case (trim(name))
    case ('bisquare', 'tukey')
      z = abs(u / max(cc, tiny(1.0_dp)))
      if (z < 1.0_dp) then
        value = 1.0_dp - (1.0_dp - z * z) ** 3
      else
        value = 1.0_dp
      end if
      if (.not. standardized) value = value * cc * cc / 6.0_dp
    case ('huber')
      z = abs(u)
      if (z <= cc) then
        raw = 0.5_dp * z * z
      else
        raw = cc * z - 0.5_dp * cc * cc
      end if
      value = raw
    case ('opt', 'optimal', 'optv0', 'mopt', 'modopt', 'moptv0')
      value = rho_shr((u / max(cc, tiny(1.0_dp))) ** 2)
      if (.not. standardized) value = 3.25_dp * cc * cc * value
    case default
      z = abs(u / max(cc, tiny(1.0_dp)))
      if (z < 1.0_dp) then
        value = 1.0_dp - (1.0_dp - z * z) ** 3
      else
        value = 1.0_dp
      end if
      if (.not. standardized) value = value * cc * cc / 6.0_dp
    end select
  end function rho_value

  elemental function rho_prime(u, family, cc, standardize) result(value)
    real(dp), intent(in) :: u, cc
    character(len=*), intent(in) :: family
    logical, intent(in), optional :: standardize
    real(dp) :: value, z, scale
    logical :: standardized
    character(len=len(family)) :: name
    name = lower_string(trim(adjustl(family)))
    standardized = .false.
    if (present(standardize)) standardized = standardize
    select case (trim(name))
    case ('bisquare', 'tukey')
      z = u / max(cc, tiny(1.0_dp))
      if (abs(z) < 1.0_dp) then
        value = u * (1.0_dp - z * z) ** 2
      else
        value = 0.0_dp
      end if
      if (standardized) value = 6.0_dp * value / (cc * cc)
    case ('huber')
      value = max(-cc, min(cc, u))
    case ('opt', 'optimal', 'optv0', 'mopt', 'modopt', 'moptv0')
      scale = max(cc, tiny(1.0_dp))
      value = 2.0_dp * u * weight_shr((u / scale) ** 2) / (scale * scale)
      if (.not. standardized) value = value * 3.25_dp * scale * scale
    case default
      z = u / max(cc, tiny(1.0_dp))
      if (abs(z) < 1.0_dp) then
        value = u * (1.0_dp - z * z) ** 2
      else
        value = 0.0_dp
      end if
      if (standardized) value = 6.0_dp * value / (cc * cc)
    end select
  end function rho_prime

  elemental function rho_second(u, family, cc, standardize) result(value)
    real(dp), intent(in) :: u, cc
    character(len=*), intent(in) :: family
    logical, intent(in), optional :: standardize
    real(dp) :: value, h
    logical :: standardized
    standardized = .false.
    if (present(standardize)) standardized = standardize
    h = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(u), abs(cc))
    value = (rho_prime(u + h, family, cc, standardized) - &
             rho_prime(u - h, family, cc, standardized)) / (2.0_dp * h)
  end function rho_second

  elemental function rho_weight(u, family, cc) result(value)
    real(dp), intent(in) :: u, cc
    character(len=*), intent(in) :: family
    real(dp) :: value
    if (abs(u) <= sqrt(epsilon(1.0_dp))) then
      value = max(rho_second(0.0_dp, family, cc), 0.0_dp)
    else
      value = rho_prime(u, family, cc) / u
    end if
    value = max(value, 0.0_dp)
  end function rho_weight

  elemental function rho_shr(d) result(value)
    real(dp), intent(in) :: d
    real(dp) :: value
    real(dp), parameter :: g1 = -1.944_dp, g2 = 1.728_dp
    real(dp), parameter :: g3 = -0.312_dp, g4 = 0.016_dp
    if (d < 4.0_dp) then
      value = 0.5_dp * max(d, 0.0_dp)
    else if (d <= 9.0_dp) then
      value = (g4 / 8.0_dp) * d ** 4 + (g3 / 6.0_dp) * d ** 3 + &
              (g2 / 4.0_dp) * d ** 2 + (g1 / 2.0_dp) * d + 1.792_dp
    else
      value = 3.25_dp
    end if
    value = value / 3.25_dp
  end function rho_shr

  elemental function weight_shr(d) result(value)
    real(dp), intent(in) :: d
    real(dp) :: value
    real(dp), parameter :: g1 = -1.944_dp, g2 = 1.728_dp
    real(dp), parameter :: g3 = -0.312_dp, g4 = 0.016_dp
    if (d < 4.0_dp) then
      value = 0.5_dp / 3.25_dp
    else if (d <= 9.0_dp) then
      value = ((g4 / 2.0_dp) * d ** 3 + (g3 / 2.0_dp) * d ** 2 + &
               (g2 / 2.0_dp) * d + g1 / 2.0_dp) / 3.25_dp
    else
      value = 0.0_dp
    end if
    value = max(value, 0.0_dp)
  end function weight_shr

  function gaussian_efficiency(family, cc) result(efficiency)
    character(len=*), intent(in) :: family
    real(dp), intent(in) :: cc
    real(dp) :: efficiency, a, b, x, h, psi, psip
    integer, parameter :: n = 4000
    integer :: i
    a = 0.0_dp
    b = 0.0_dp
    h = 10.0_dp / real(n, dp)
    do i = 0, n
      x = h * real(i, dp)
      psi = rho_prime(x, family, cc)
      psip = rho_second(x, family, cc)
      if (i == 0 .or. i == n) then
        a = a + psi * psi * normal_density(x)
        b = b + psip * normal_density(x)
      else if (mod(i, 2) == 1) then
        a = a + 4.0_dp * psi * psi * normal_density(x)
        b = b + 4.0_dp * psip * normal_density(x)
      else
        a = a + 2.0_dp * psi * psi * normal_density(x)
        b = b + 2.0_dp * psip * normal_density(x)
      end if
    end do
    a = 2.0_dp * h * a / 3.0_dp
    b = 2.0_dp * h * b / 3.0_dp
    if (a <= tiny(1.0_dp)) then
      efficiency = 0.0_dp
    else
      efficiency = max(0.0_dp, min(1.0_dp, b * b / a))
    end if
  end function gaussian_efficiency

  function tuning_for_efficiency(efficiency, family) result(cc)
    real(dp), intent(in) :: efficiency
    character(len=*), intent(in) :: family
    real(dp) :: cc, lo, hi, mid, emid, target
    integer :: i
    target = max(0.05_dp, min(0.9999_dp, efficiency))
    lo = 0.05_dp
    hi = 25.0_dp
    do i = 1, 80
      mid = 0.5_dp * (lo + hi)
      emid = gaussian_efficiency(family, mid)
      if (emid < target) then
        lo = mid
      else
        hi = mid
      end if
    end do
    cc = 0.5_dp * (lo + hi)
  end function tuning_for_efficiency

  function tuning_huber(efficiency) result(cc)
    real(dp), intent(in) :: efficiency
    real(dp) :: cc
    cc = tuning_for_efficiency(efficiency, 'huber')
  end function tuning_huber

  function tuning_bisquare(efficiency) result(cc)
    real(dp), intent(in) :: efficiency
    real(dp) :: cc
    cc = tuning_for_efficiency(efficiency, 'bisquare')
  end function tuning_bisquare

  function tuning_opt(efficiency) result(cc)
    real(dp), intent(in) :: efficiency
    real(dp) :: cc
    cc = tuning_for_efficiency(efficiency, 'opt')
  end function tuning_opt

  function tuning_mopt(efficiency) result(cc)
    real(dp), intent(in) :: efficiency
    real(dp) :: cc
    cc = tuning_for_efficiency(efficiency, 'mopt')
  end function tuning_mopt

  function scale_m(residuals, delta, family, tuning, max_iter, tol) result(scale)
    real(dp), intent(in) :: residuals(:)
    real(dp), intent(in), optional :: delta, tuning, tol
    character(len=*), intent(in), optional :: family
    integer, intent(in), optional :: max_iter
    real(dp) :: scale, d, c, tolerance, next_scale, error_value
    character(len=16) :: fam
    integer :: iteration, iterations
    d = 0.5_dp
    if (present(delta)) d = delta
    fam = 'bisquare'
    if (present(family)) fam = family
    c = 1.54764_dp
    if (present(tuning)) c = tuning
    iterations = 100
    if (present(max_iter)) iterations = max(1, max_iter)
    tolerance = 1.0e-6_dp
    if (present(tol)) tolerance = max(tol, epsilon(1.0_dp))
    if (size(residuals) == 0) then
      scale = 0.0_dp
      return
    end if
    scale = median_absolute(residuals) / 0.6745_dp
    if (scale <= epsilon(1.0_dp)) then
      scale = 0.0_dp
      return
    end if
    do iteration = 1, iterations
      next_scale = scale * sqrt(max(mean_rho(residuals / scale, fam, c) / d, 0.0_dp))
      error_value = abs(next_scale - scale) / max(scale, tiny(1.0_dp))
      scale = next_scale
      if (error_value <= tolerance) exit
    end do
  end function scale_m

  function mean_rho(x, family, cc) result(value)
    real(dp), intent(in) :: x(:), cc
    character(len=*), intent(in) :: family
    real(dp) :: value
    integer :: i
    value = 0.0_dp
    do i = 1, size(x)
      value = value + rho_value(x(i), family, cc, .true.)
    end do
    value = value / real(max(1, size(x)), dp)
  end function mean_rho

  function expected_rho(scale, family, cc) result(value)
    real(dp), intent(in) :: scale, cc
    character(len=*), intent(in) :: family
    real(dp) :: value, x, h, term
    integer, parameter :: n = 4000
    integer :: i
    h = 10.0_dp / real(n, dp)
    value = 0.0_dp
    do i = 0, n
      x = h * real(i, dp)
      term = rho_value(x / max(scale, tiny(1.0_dp)), family, cc, .true.) * normal_density(x)
      if (i == 0 .or. i == n) then
        value = value + term
      else if (mod(i, 2) == 1) then
        value = value + 4.0_dp * term
      else
        value = value + 2.0_dp * term
      end if
    end do
    value = 2.0_dp * h * value / 3.0_dp
  end function expected_rho

  function transformed_r_squared(r2, family, cc) result(value)
    real(dp), intent(in) :: r2, cc
    character(len=*), intent(in) :: family
    real(dp) :: value, a, b
    a = expected_rho(1.0_dp, family, cc)
    b = expected_rho(sqrt(max(1.0_dp - r2, 1.0e-12_dp)), family, cc)
    value = (b - a) / max(b * (1.0_dp - a), tiny(1.0_dp))
  end function transformed_r_squared

  function inverse_robust_r_squared(rr2, family, cc) result(r2)
    real(dp), intent(in) :: rr2, cc
    character(len=*), intent(in) :: family
    real(dp) :: r2, lo, hi, mid, value
    integer :: i
    if (rr2 >= 0.99_dp) then
      r2 = 1.0_dp
      return
    end if
    if (rr2 <= transformed_r_squared(1.0e-5_dp, family, cc)) then
      r2 = 0.0_dp
      return
    end if
    lo = 0.0_dp
    hi = 0.999999_dp
    do i = 1, 80
      mid = 0.5_dp * (lo + hi)
      value = transformed_r_squared(mid, family, cc)
      if (value < rr2) then
        lo = mid
      else
        hi = mid
      end if
    end do
    r2 = 0.5_dp * (lo + hi)
  end function inverse_robust_r_squared
end module robstattm_psi
