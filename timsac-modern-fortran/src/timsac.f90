! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

module timsac
  use timsac_kinds, only: dp
  use timsac_raw, only: autcorf, mulcorf, auspecf, wnoisef
  implicit none
  private

  integer, parameter, public :: window_hanning = 1
  integer, parameter, public :: window_akaike = 2

  type, public :: autocorrelation_result
    real(dp) :: mean = 0.0_dp
    real(dp), allocatable :: covariance(:)
    real(dp), allocatable :: correlation(:)
  end type autocorrelation_result

  type, public :: multivariate_correlation_result
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: covariance(:,:,:)
    real(dp), allocatable :: correlation(:,:,:)
  end type multivariate_correlation_result

  type, public :: spectrum_result
    real(dp), allocatable :: frequency(:)
    real(dp), allocatable :: spectrum(:)
    real(dp), allocatable :: statistic(:)
  end type spectrum_result

  public :: dp
  public :: autocorrelation
  public :: multivariate_correlation
  public :: power_spectrum
  public :: matrix_filter
  public :: white_noise

contains

  function autocorrelation(x, max_lag) result(out)
    !! Compute biased autocovariances and autocorrelations.
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: max_lag
    type(autocorrelation_result) :: out

    integer :: lag, lag1, n
    real(dp), allocatable :: work(:)

    n = size(x)
    if (n < 1) error stop "autocorrelation: x must not be empty"

    lag = default_lag(n)
    if (present(max_lag)) lag = max_lag
    if (lag < 0 .or. lag >= n) then
      error stop "autocorrelation: max_lag must be between 0 and size(x)-1"
    end if

    lag1 = lag + 1
    allocate(work(n))
    work = x
    allocate(out%covariance(0:lag), out%correlation(0:lag))

    call autcorf(work, n, out%covariance, out%correlation, lag1, out%mean)
  end function autocorrelation


  function multivariate_correlation(x, max_lag) result(out)
    !! Compute biased auto- and cross-covariances for a multivariate series.
    real(dp), intent(in) :: x(:,:)
    integer, intent(in), optional :: max_lag
    type(multivariate_correlation_result) :: out

    integer :: d, lag, lag1, n
    real(dp), allocatable :: work(:,:)

    n = size(x, 1)
    d = size(x, 2)
    if (n < 1 .or. d < 1) then
      error stop "multivariate_correlation: x must have positive dimensions"
    end if

    lag = default_lag(n)
    if (present(max_lag)) lag = max_lag
    if (lag < 0 .or. lag >= n) then
      error stop "multivariate_correlation: max_lag must be between 0 and size(x,1)-1"
    end if

    lag1 = lag + 1
    allocate(work(n, d))
    work = x
    allocate(out%mean(d))
    allocate(out%covariance(0:lag, d, d))
    allocate(out%correlation(0:lag, d, d))

    call mulcorf(work, n, d, lag1, out%mean, out%covariance, out%correlation)
  end function multivariate_correlation


  function power_spectrum(x, max_lag, window) result(out)
    !! Estimate a univariate spectrum using a TIMSAC Blackman-Tukey window.
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: max_lag
    integer, intent(in), optional :: window
    type(spectrum_result) :: out

    integer :: i, lag, lag1, n, selected_window
    real(dp), allocatable :: p_hanning(:), p_akaike(:)
    type(autocorrelation_result) :: ac

    n = size(x)
    if (n < 2) error stop "power_spectrum: at least two observations are required"

    lag = default_lag(n)
    if (present(max_lag)) lag = max_lag
    if (lag < 1 .or. lag >= n) then
      error stop "power_spectrum: max_lag must be between 1 and size(x)-1"
    end if

    selected_window = window_akaike
    if (present(window)) selected_window = window
    if (selected_window /= window_hanning .and. selected_window /= window_akaike) then
      error stop "power_spectrum: invalid window"
    end if

    ac = autocorrelation(x, lag)
    lag1 = lag + 1
    allocate(p_hanning(0:lag), p_akaike(0:lag))
    allocate(out%frequency(0:lag), out%spectrum(0:lag), out%statistic(0:lag))

    call auspecf(n, lag1, ac%covariance, p_hanning, p_akaike, out%statistic)
    if (selected_window == window_hanning) then
      out%spectrum = p_hanning
    else
      out%spectrum = p_akaike
    end if

    do i = 0, lag
      out%frequency(i) = real(i, dp) / real(2 * lag, dp)
    end do
  end function power_spectrum


  function matrix_filter(x, coefficients, recursive, initial) result(y)
    !! Apply the multivariate filter used by the R package's mfilter().
    !!
    !! coefficients(:,:,j) is the j-th matrix lag. With recursive=.false.,
    !! y(t) = x(t) - sum_j coefficients(:,:,j) * x(t-j). With
    !! recursive=.true., y(t) = x(t) + sum_j coefficients(:,:,j) * y(t-j).
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in) :: coefficients(:,:,:)
    logical, intent(in), optional :: recursive
    real(dp), intent(in), optional :: initial(:,:)
    real(dp), allocatable :: y(:,:)

    integer :: d, j, n, p, t
    logical :: use_recursive
    real(dp), allocatable :: history(:,:)

    n = size(x, 1)
    d = size(x, 2)
    p = size(coefficients, 3)

    if (n < 1 .or. d < 1 .or. p < 1) then
      error stop "matrix_filter: x and coefficients must have positive dimensions"
    end if
    if (size(coefficients, 1) /= d .or. size(coefficients, 2) /= d) then
      error stop "matrix_filter: coefficient matrices must be d by d"
    end if
    if (present(initial)) then
      if (size(initial, 1) /= p .or. size(initial, 2) /= d) then
        error stop "matrix_filter: initial must have shape [p,d]"
      end if
    end if

    use_recursive = .false.
    if (present(recursive)) use_recursive = recursive

    allocate(y(n, d))
    allocate(history(p + n, d))
    history = 0.0_dp
    if (present(initial)) history(1:p, :) = initial

    if (use_recursive) then
      do t = 1, n
        history(p + t, :) = x(t, :)
        do j = 1, p
          history(p + t, :) = history(p + t, :) + &
            matmul(coefficients(:, :, j), history(p + t - j, :))
        end do
        y(t, :) = history(p + t, :)
      end do
    else
      do t = 1, n
        y(t, :) = x(t, :)
        do j = 1, p
          y(t, :) = y(t, :) - &
            matmul(coefficients(:, :, j), history(p + t - j, :))
        end do
        history(p + t, :) = x(t, :)
      end do
    end if
  end function matrix_filter


  subroutine white_noise(covariance, n, noise)
    !! Generate approximately Gaussian vector white noise.
    real(dp), intent(in) :: covariance(:,:)
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: noise(:,:)

    integer :: d
    real(dp), allocatable :: covariance_work(:,:)

    d = size(covariance, 1)
    if (d < 1 .or. size(covariance, 2) /= d) then
      error stop "white_noise: covariance must be a nonempty square matrix"
    end if
    if (n < 1) error stop "white_noise: n must be positive"

    allocate(covariance_work(d, d), noise(d, n))
    covariance_work = covariance
    call wnoisef(n, d, covariance_work, noise)
  end subroutine white_noise


  integer function default_lag(n) result(lag)
    integer, intent(in) :: n

    lag = int(2.0_dp * sqrt(real(n, dp)))
    lag = min(lag, n - 1)
  end function default_lag

end module timsac
