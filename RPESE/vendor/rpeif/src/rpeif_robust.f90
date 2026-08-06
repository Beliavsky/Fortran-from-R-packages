! SPDX-License-Identifier: GPL-3.0-or-later
module rpeif_robust
  use rpeif_kinds, only : dp
  use rpeif_stats, only : median_value, mad_value, mean_value, normal_quantile, lower_string
  use robstattm_psi, only : rho_prime, rho_second, rho_weight, tuning_for_efficiency, scale_m
  implicit none
  private
  public :: robust_location_scale, robust_mean_influence, robust_clean
contains
  subroutine robust_location_scale(x, location, scale, family, efficiency, converged)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: location, scale
    character(len=*), intent(in), optional :: family
    real(dp), intent(in), optional :: efficiency
    logical, intent(out), optional :: converged
    character(len=16) :: fam
    real(dp) :: eff, tuning, old_location, new_location, initial_scale
    real(dp), allocatable :: residuals(:), weights(:)
    integer :: i, iter
    logical :: ok

    fam = 'mopt'
    if (present(family)) fam = lower_string(trim(adjustl(family)))
    eff = 0.95_dp
    if (present(efficiency)) eff = efficiency
    ok = .false.

    if (size(x) == 0) then
      location = 0.0_dp
      scale = 0.0_dp
      if (present(converged)) converged = .false.
      return
    end if

    tuning = tuning_for_efficiency(eff, fam)
    old_location = median_value(x)
    initial_scale = mad_value(x)
    if (initial_scale <= 1.0e-14_dp) then
      location = old_location
      scale = 0.0_dp
      if (present(converged)) converged = .true.
      return
    end if

    allocate(residuals(size(x)), weights(size(x)))
    do iter = 1, 100
      residuals = (x - old_location) / initial_scale
      do i = 1, size(x)
        weights(i) = rho_weight(residuals(i), fam, tuning)
      end do
      if (sum(weights) <= tiny(1.0_dp)) exit
      new_location = sum(weights * x) / sum(weights)
      if (abs(new_location - old_location) <= 1.0e-7_dp * initial_scale) then
        old_location = new_location
        ok = .true.
        exit
      end if
      old_location = new_location
    end do

    location = old_location
    scale = scale_m(x - location, 0.5_dp, fam, tuning, 100, 1.0e-7_dp)
    if (scale <= 0.0_dp) scale = initial_scale
    if (present(converged)) converged = ok
  end subroutine robust_location_scale

  subroutine robust_mean_influence(x_eval, returns, values, family, efficiency, status)
    real(dp), intent(in) :: x_eval(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    character(len=*), intent(in), optional :: family
    real(dp), intent(in), optional :: efficiency
    integer, intent(out), optional :: status
    character(len=16) :: fam
    real(dp) :: eff, tuning, location, scale, denominator
    real(dp), allocatable :: standardized(:)
    integer :: i
    logical :: converged

    fam = 'mopt'
    if (present(family)) fam = lower_string(trim(adjustl(family)))
    eff = 0.95_dp
    if (present(efficiency)) eff = efficiency
    allocate(values(size(x_eval)))
    values = 0.0_dp
    if (size(returns) == 0) then
      if (present(status)) status = 1
      return
    end if

    call robust_location_scale(returns, location, scale, fam, eff, converged)
    scale = mad_value(returns)
    if (scale <= 1.0e-14_dp) then
      if (present(status)) status = 2
      return
    end if
    tuning = tuning_for_efficiency(eff, fam)
    allocate(standardized(size(returns)))
    standardized = (returns - location) / scale
    denominator = 0.0_dp
    do i = 1, size(returns)
      denominator = denominator + rho_second(standardized(i), fam, tuning)
    end do
    denominator = denominator / real(size(returns), dp)
    if (abs(denominator) <= 1.0e-14_dp) then
      if (present(status)) status = 2
      return
    end if
    do i = 1, size(x_eval)
      values(i) = scale * rho_prime((x_eval(i) - location) / scale, fam, tuning) / denominator
    end do
    if (present(status)) status = merge(0, 2, converged)
  end subroutine robust_mean_influence

  subroutine robust_clean(x, cleaned, efficiency, family, status)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: cleaned(:)
    real(dp), intent(in), optional :: efficiency
    character(len=*), intent(in), optional :: family
    integer, intent(out), optional :: status
    real(dp) :: eff, location, scale, cutoff
    character(len=16) :: fam
    logical :: converged

    eff = 0.95_dp
    if (present(efficiency)) eff = efficiency
    fam = 'mopt'
    if (present(family)) fam = lower_string(trim(adjustl(family)))
    call robust_location_scale(x, location, scale, fam, eff, converged)
    allocate(cleaned(size(x)))
    if (scale <= 0.0_dp) then
      cleaned = x
      if (present(status)) status = 2
      return
    end if

    if (abs(eff - 0.95_dp) <= 1.0e-12_dp) then
      cutoff = 3.0_dp
    else if (abs(eff - 0.99_dp) <= 1.0e-12_dp) then
      cutoff = 3.568_dp
    else if (abs(eff - 0.999_dp) <= 1.0e-12_dp) then
      cutoff = 4.21_dp
    else
      cutoff = normal_quantile(0.5_dp * (1.0_dp + max(0.0_dp, min(0.999999_dp, eff))))
      cutoff = max(cutoff, 1.0_dp)
    end if
    cleaned = max(location - cutoff * scale, min(location + cutoff * scale, x))
    if (present(status)) status = merge(0, 2, converged)
  end subroutine robust_clean
end module rpeif_robust
