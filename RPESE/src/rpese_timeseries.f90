! SPDX-License-Identifier: GPL-3.0-or-later
module rpese_timeseries
  use rpese_kinds, only : dp, pi
  use rpese_types, only : periodogram_result, rpese_options, rpese_success, &
    rpese_invalid_argument, rpese_numerical_failure, frequency_all, &
    frequency_decimate, frequency_truncate
  use rpeif_stats, only : mean_value, sample_sd
  implicit none
  private
  public :: fit_periodogram, lag1_correlation, fit_ar1, polynomial_design
contains
  subroutine fit_periodogram(data, result, options, twosided)
    real(dp), intent(in) :: data(:)
    type(periodogram_result), intent(out) :: result
    type(rpese_options), intent(in), optional :: options
    logical, intent(in), optional :: twosided
    type(rpese_options) :: opts
    logical :: two
    real(dp), allocatable :: full_frequency(:), full_spectrum(:)
    real(dp) :: angle, real_part, imag_part
    integer :: n, m, k, t, nkeep, stride, i, j

    opts = rpese_options()
    if (present(options)) opts = options
    two = .false.
    if (present(twosided)) two = twosided
    n = size(data)
    m = n / 2 - 1
    if (n < 4 .or. m < 1) then
      result%status = rpese_invalid_argument
      result%message = 'At least four observations are required for a periodogram.'
      allocate(result%frequency(0), result%spectrum(0))
      return
    end if
    if (opts%frequency_fraction <= 0.0_dp .or. opts%frequency_fraction > 1.0_dp .or. &
        opts%keep_fraction <= 0.0_dp .or. opts%keep_fraction > 1.0_dp) then
      result%status = rpese_invalid_argument
      result%message = 'Frequency and keep fractions must lie in (0, 1].'
      allocate(result%frequency(0), result%spectrum(0))
      return
    end if

    allocate(full_frequency(m), full_spectrum(m))
    do k = 1, m
      real_part = 0.0_dp
      imag_part = 0.0_dp
      do t = 1, n
        angle = 2.0_dp * pi * real(k * (t - 1), dp) / real(n, dp)
        real_part = real_part + data(t) * cos(angle)
        imag_part = imag_part - data(t) * sin(angle)
      end do
      full_spectrum(k) = max(1.0e-5_dp, (real_part * real_part + imag_part * imag_part) / real(n, dp))
      full_frequency(k) = real(k, dp) / real(n, dp)
    end do

    select case (opts%frequency_mode)
    case (frequency_all)
      nkeep = max(1, min(m, int(floor(real(m, dp) * opts%keep_fraction))))
      allocate(result%frequency(nkeep), result%spectrum(nkeep))
      result%frequency = full_frequency(1:nkeep)
      result%spectrum = full_spectrum(1:nkeep)
    case (frequency_decimate)
      stride = max(1, nint(1.0_dp / opts%frequency_fraction))
      nkeep = 1 + (m - 1) / stride
      allocate(result%frequency(nkeep), result%spectrum(nkeep))
      j = 0
      do i = 1, m, stride
        j = j + 1
        result%frequency(j) = full_frequency(i)
        result%spectrum(j) = full_spectrum(i)
      end do
    case (frequency_truncate)
      nkeep = max(1, min(m, int(floor(real(m, dp) * opts%frequency_fraction))))
      allocate(result%frequency(nkeep), result%spectrum(nkeep))
      result%frequency = full_frequency(1:nkeep)
      result%spectrum = full_spectrum(1:nkeep)
    case default
      result%status = rpese_invalid_argument
      result%message = 'Unknown frequency inclusion mode.'
      allocate(result%frequency(0), result%spectrum(0))
      return
    end select

    if (two) call make_two_sided(result)
    result%status = rpese_success
    result%message = 'completed'
  end subroutine fit_periodogram

  subroutine make_two_sided(result)
    type(periodogram_result), intent(inout) :: result
    real(dp), allocatable :: frequency(:), spectrum(:)
    integer :: n, i
    n = size(result%frequency)
    allocate(frequency(2 * n), spectrum(2 * n))
    do i = 1, n
      frequency(i) = -result%frequency(n - i + 1)
      spectrum(i) = result%spectrum(n - i + 1)
      frequency(n + i) = result%frequency(i)
      spectrum(n + i) = result%spectrum(i)
    end do
    call move_alloc(frequency, result%frequency)
    call move_alloc(spectrum, result%spectrum)
  end subroutine make_two_sided

  real(dp) function lag1_correlation(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: mean_left, mean_right, denominator_left, denominator_right
    if (size(x) < 3) then
      value = 0.0_dp
      return
    end if
    mean_left = sum(x(1:size(x)-1)) / real(size(x) - 1, dp)
    mean_right = sum(x(2:size(x))) / real(size(x) - 1, dp)
    denominator_left = sum((x(1:size(x)-1) - mean_left) ** 2)
    denominator_right = sum((x(2:size(x)) - mean_right) ** 2)
    if (denominator_left <= tiny(1.0_dp) .or. denominator_right <= tiny(1.0_dp)) then
      value = 0.0_dp
    else
      value = sum((x(1:size(x)-1) - mean_left) * (x(2:size(x)) - mean_right)) / &
        sqrt(denominator_left * denominator_right)
    end if
  end function lag1_correlation

  subroutine fit_ar1(x, coefficient, intercept, residuals, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: coefficient, intercept
    real(dp), allocatable, intent(out), optional :: residuals(:)
    integer, intent(out), optional :: status
    real(dp) :: xbar, ybar, denominator, unconditional_mean
    integer :: n

    n = size(x)
    coefficient = 0.0_dp
    intercept = 0.0_dp
    if (n < 3) then
      if (present(residuals)) then
        allocate(residuals(n))
        residuals = x - mean_value(x)
      end if
      if (present(status)) status = rpese_invalid_argument
      return
    end if
    xbar = sum(x(1:n-1)) / real(n - 1, dp)
    ybar = sum(x(2:n)) / real(n - 1, dp)
    denominator = sum((x(1:n-1) - xbar) ** 2)
    if (denominator <= tiny(1.0_dp)) then
      if (present(residuals)) then
        allocate(residuals(n))
        residuals = x - mean_value(x)
      end if
      if (present(status)) status = rpese_numerical_failure
      return
    end if
    coefficient = sum((x(1:n-1) - xbar) * (x(2:n) - ybar)) / denominator
    coefficient = max(-0.999_dp, min(0.999_dp, coefficient))
    intercept = ybar - coefficient * xbar
    if (present(residuals)) then
      allocate(residuals(n))
      if (abs(1.0_dp - coefficient) > 1.0e-8_dp) then
        unconditional_mean = intercept / (1.0_dp - coefficient)
      else
        unconditional_mean = mean_value(x)
      end if
      residuals(1) = x(1) - unconditional_mean
      residuals(2:n) = x(2:n) - intercept - coefficient * x(1:n-1)
    end if
    if (present(status)) status = rpese_success
  end subroutine fit_ar1

  subroutine polynomial_design(frequency, degree, design, standardize, means, scales, status)
    real(dp), intent(in) :: frequency(:)
    integer, intent(in) :: degree
    real(dp), allocatable, intent(out) :: design(:, :)
    logical, intent(in), optional :: standardize
    real(dp), allocatable, intent(out), optional :: means(:), scales(:)
    integer, intent(out), optional :: status
    logical :: do_standardize
    real(dp), allocatable :: local_means(:), local_scales(:)
    integer :: j

    do_standardize = .false.
    if (present(standardize)) do_standardize = standardize
    if (size(frequency) < 1 .or. degree < 0) then
      allocate(design(0, 0), local_means(0), local_scales(0))
      if (present(means)) means = local_means
      if (present(scales)) scales = local_scales
      if (present(status)) status = rpese_invalid_argument
      return
    end if
    allocate(design(size(frequency), degree + 1))
    design(:, 1) = 1.0_dp
    do j = 1, degree
      design(:, j + 1) = frequency ** j
    end do
    allocate(local_means(degree), local_scales(degree))
    do j = 1, degree
      local_means(j) = mean_value(design(:, j + 1))
      local_scales(j) = sample_sd(design(:, j + 1))
      if (local_scales(j) <= tiny(1.0_dp)) local_scales(j) = 1.0_dp
      if (do_standardize) design(:, j + 1) = &
        (design(:, j + 1) - local_means(j)) / local_scales(j)
    end do
    if (present(means)) means = local_means
    if (present(scales)) scales = local_scales
    if (present(status)) status = rpese_success
  end subroutine polynomial_design
end module rpese_timeseries
