! SPDX-License-Identifier: GPL-2.0-only
module tsmarch_fft
  use ghyp_kinds, only : dp
  use tsmarch_types, only : fft_distribution, tsm_success, tsm_invalid_argument
  implicit none
  private
  public :: discrete_convolution, fft_convolution, dfft, pfft, qfft

  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  function discrete_convolution(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp), allocatable :: c(:)
    integer :: i, j
    allocate(c(size(a) + size(b) - 1))
    c = 0.0_dp
    do j = 1, size(b)
      do i = 1, size(a)
        c(i + j - 1) = c(i + j - 1) + a(i) * b(j)
      end do
    end do
  end function discrete_convolution

  function fft_convolution(component_density, lower, step) result(out)
    real(dp), intent(in) :: component_density(:, :)
    real(dp), intent(in) :: lower, step
    type(fft_distribution) :: out
    complex(dp), allocatable :: spectrum(:), transform(:)
    real(dp), allocatable :: density(:), work(:)
    integer :: n, m, i, j, k
    real(dp) :: total
    n = size(component_density, 1)
    m = size(component_density, 2)
    if (n < 2 .or. m < 1 .or. step <= 0.0_dp .or. any(component_density < 0.0_dp)) then
      out%status = tsm_invalid_argument
      out%message = 'invalid component densities or grid spacing'
      return
    end if
    ! A zero-padded direct DFT is used for deterministic portability.
    k = 1
    do while (k < m * (n - 1) + 1)
      k = 2 * k
    end do
    allocate(spectrum(k), transform(k), work(k))
    spectrum = cmplx(1.0_dp, 0.0_dp, dp)
    do j = 1, m
      work = 0.0_dp
      work(1:n) = component_density(:, j)
      transform = dft_real(work)
      spectrum = spectrum * transform
    end do
    work = inverse_dft_real(spectrum)
    allocate(density(m * (n - 1) + 1))
    density = max(work(1:size(density)), 0.0_dp)
    total = sum(density) * step
    if (total > tiny(1.0_dp)) density = density / total
    allocate(out%grid(size(density)), out%density(size(density)), out%cdf(size(density)))
    out%grid = [(real(i - 1, dp) * step + real(m, dp) * lower, i = 1, size(density))]
    out%density = density
    out%cdf(1) = 0.5_dp * density(1) * step
    do i = 2, size(density)
      out%cdf(i) = out%cdf(i - 1) + 0.5_dp * step * (density(i - 1) + density(i))
    end do
    if (out%cdf(size(out%cdf)) > 0.0_dp) out%cdf = out%cdf / out%cdf(size(out%cdf))
    out%status = tsm_success
    out%message = 'ok'
  end function fft_convolution

  function dfft(object, x) result(value)
    type(fft_distribution), intent(in) :: object
    real(dp), intent(in) :: x
    real(dp) :: value
    value = interpolate(object%grid, object%density, x, 0.0_dp, 0.0_dp)
  end function dfft

  function pfft(object, x, lower_tail) result(value)
    type(fft_distribution), intent(in) :: object
    real(dp), intent(in) :: x
    logical, intent(in), optional :: lower_tail
    real(dp) :: value
    logical :: lower
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    value = interpolate(object%grid, object%cdf, x, 0.0_dp, 1.0_dp)
    if (.not. lower) value = 1.0_dp - value
  end function pfft

  function qfft(object, probability, lower_tail) result(value)
    type(fft_distribution), intent(in) :: object
    real(dp), intent(in) :: probability
    logical, intent(in), optional :: lower_tail
    real(dp) :: value, p, frac
    logical :: lower
    integer :: i
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    p = probability
    if (.not. lower) p = 1.0_dp - p
    p = min(max(p, 0.0_dp), 1.0_dp)
    if (p <= object%cdf(1)) then
      value = object%grid(1)
      return
    end if
    do i = 2, size(object%cdf)
      if (object%cdf(i) >= p) then
        frac = (p - object%cdf(i - 1)) / max(object%cdf(i) - object%cdf(i - 1), tiny(1.0_dp))
        value = object%grid(i - 1) + frac * (object%grid(i) - object%grid(i - 1))
        return
      end if
    end do
    value = object%grid(size(object%grid))
  end function qfft

  function dft_real(x) result(out)
    real(dp), intent(in) :: x(:)
    complex(dp), allocatable :: out(:)
    complex(dp) :: phase
    integer :: n, k, j
    n = size(x)
    allocate(out(n))
    do k = 1, n
      out(k) = cmplx(0.0_dp, 0.0_dp, dp)
      do j = 1, n
        phase = exp(cmplx(0.0_dp, -2.0_dp * pi * real((k - 1) * (j - 1), dp) / real(n, dp), dp))
        out(k) = out(k) + x(j) * phase
      end do
    end do
  end function dft_real

  function inverse_dft_real(x) result(out)
    complex(dp), intent(in) :: x(:)
    real(dp), allocatable :: out(:)
    complex(dp) :: phase, value
    integer :: n, k, j
    n = size(x)
    allocate(out(n))
    do j = 1, n
      value = cmplx(0.0_dp, 0.0_dp, dp)
      do k = 1, n
        phase = exp(cmplx(0.0_dp, 2.0_dp * pi * real((k - 1) * (j - 1), dp) / real(n, dp), dp))
        value = value + x(k) * phase
      end do
      out(j) = real(value, dp) / real(n, dp)
    end do
  end function inverse_dft_real

  function interpolate(x, y, z, left_value, right_value) result(value)
    real(dp), intent(in) :: x(:), y(:), z, left_value, right_value
    real(dp) :: value, frac
    integer :: i
    if (z <= x(1)) then
      value = left_value
      return
    end if
    if (z >= x(size(x))) then
      value = right_value
      return
    end if
    do i = 2, size(x)
      if (x(i) >= z) then
        frac = (z - x(i - 1)) / (x(i) - x(i - 1))
        value = y(i - 1) + frac * (y(i) - y(i - 1))
        return
      end if
    end do
    value = right_value
  end function interpolate

end module tsmarch_fft
