! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_bias
  use multiatsm_kinds, only : dp
  use multiatsm_linalg, only : spectral_radius
  use multiatsm_random, only : set_random_seed, random_integer
  use multiatsm_types, only : var_model
  use multiatsm_var, only : fit_var
  use multiatsm_bootstrap, only : simulate_var
  implicit none
  private

  public :: bias_correct_var, shrink_transition

contains

  subroutine bias_correct_var(data, n_bootstrap, burn_iterations, averaging_iterations, gamma, &
      corrected_phi, corrected_sigma, info, seed)
    real(dp), intent(in) :: data(:, :)
    integer, intent(in) :: n_bootstrap, burn_iterations, averaging_iterations
    real(dp), intent(in) :: gamma
    real(dp), allocatable, intent(out) :: corrected_phi(:, :), corrected_sigma(:, :)
    integer, intent(out) :: info
    integer, intent(in), optional :: seed
    type(var_model) :: original, fitted
    real(dp), allocatable :: current(:, :), mean_estimate(:, :), sum_after_burn(:, :)
    real(dp), allocatable :: residual_draw(:, :), simulated(:, :), mean_data(:), intercept(:)
    real(dp), allocatable :: sigma_sum(:, :)
    integer :: k, t, iteration, b, total_iterations, local_seed, j

    if (n_bootstrap < 1 .or. burn_iterations < 0 .or. averaging_iterations < 1 .or. &
        gamma <= 0.0_dp .or. gamma > 1.0_dp) then
      allocate(corrected_phi(0, 0), corrected_sigma(0, 0))
      info = -1
      return
    end if
    call fit_var(data, original, info)
    if (info /= 0) return
    k = size(data, 1)
    t = size(data, 2)
    allocate(current(k, k), mean_estimate(k, k), sum_after_burn(k, k), sigma_sum(k, k))
    allocate(residual_draw(k, t - 1), mean_data(k), intercept(k))
    current = original%phi
    sum_after_burn = 0.0_dp
    sigma_sum = 0.0_dp
    mean_data = sum(data, dim=2) / real(t, dp)
    local_seed = 24680
    if (present(seed)) local_seed = seed
    total_iterations = burn_iterations + averaging_iterations

    do iteration = 1, total_iterations
      mean_estimate = 0.0_dp
      sigma_sum = 0.0_dp
      intercept = mean_data - matmul(current, mean_data)
      do b = 1, n_bootstrap
        call set_random_seed(local_seed + 100003 * iteration + 7919 * b)
        do j = 1, t - 1
          residual_draw(:, j) = original%residuals(:, random_integer(t - 1))
        end do
        call simulate_var(intercept, current, residual_draw, data(:, random_integer(t)), simulated, info)
        if (info /= 0) return
        call fit_var(simulated, fitted, info)
        if (info /= 0) return
        mean_estimate = mean_estimate + fitted%phi
        sigma_sum = sigma_sum + fitted%sigma
      end do
      mean_estimate = mean_estimate / real(n_bootstrap, dp)
      current = current + gamma * (original%phi - mean_estimate)
      if (iteration > burn_iterations) sum_after_burn = sum_after_burn + current
    end do
    allocate(corrected_phi(k, k), corrected_sigma(k, k))
    corrected_phi = sum_after_burn / real(averaging_iterations, dp)
    corrected_sigma = sigma_sum / real(n_bootstrap, dp)
    info = 0
  end subroutine bias_correct_var

  subroutine shrink_transition(corrected, uncorrected, maximum_radius, shrunk, info, tolerance)
    real(dp), intent(in) :: corrected(:, :), uncorrected(:, :), maximum_radius
    real(dp), allocatable, intent(out) :: shrunk(:, :)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: tolerance
    real(dp) :: radius_corrected, radius_uncorrected, radius_trial, lo, hi, mid, tol
    real(dp), allocatable :: trial(:, :)
    integer :: iteration

    if (size(corrected, 1) /= size(corrected, 2) .or. any(shape(corrected) /= shape(uncorrected)) .or. &
        maximum_radius <= 0.0_dp) then
      allocate(shrunk(0, 0))
      info = -1
      return
    end if
    tol = 1.0e-8_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    radius_corrected = spectral_radius(corrected, info)
    if (info /= 0) return
    radius_uncorrected = spectral_radius(uncorrected, info)
    if (info /= 0) return
    allocate(shrunk(size(corrected, 1), size(corrected, 2)))
    if (radius_corrected <= maximum_radius) then
      shrunk = corrected
      info = 0
      return
    end if
    if (radius_uncorrected >= maximum_radius) then
      shrunk = uncorrected
      info = 1
      return
    end if
    allocate(trial(size(corrected, 1), size(corrected, 2)))
    lo = 0.0_dp
    hi = 1.0_dp
    do iteration = 1, 100
      mid = 0.5_dp * (lo + hi)
      trial = corrected + mid * (uncorrected - corrected)
      radius_trial = spectral_radius(trial, info)
      if (info /= 0) return
      if (radius_trial > maximum_radius) then
        lo = mid
      else
        hi = mid
      end if
      if (hi - lo <= tol) exit
    end do
    shrunk = corrected + hi * (uncorrected - corrected)
    info = 0
  end subroutine shrink_transition

end module multiatsm_bias
