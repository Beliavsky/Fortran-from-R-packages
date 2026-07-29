! SPDX-License-Identifier: GPL-3.0-or-later
module frapo_series
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use frapo_kinds, only : dp
  use frapo_types, only : frapo_ok, frapo_invalid_input
  use frapo_linalg, only : solve_linear
  implicit none
  private

  integer, parameter, public :: returns_continuous = 1
  integer, parameter, public :: returns_discrete = 2
  integer, parameter, public :: convert_cont_to_disc = 1
  integer, parameter, public :: convert_disc_to_cont = 2

  public :: cap_series, return_series, return_convert
  public :: trend_binary, trend_bilson, trend_exponential
  public :: trend_hodrick_prescott, trend_simple_moving_average
  public :: trend_weighted_moving_average
  public :: capser, returnseries, returnconvert
  public :: trdbinary, trdbilson, trdes, trdhp, trdsma, trdwma

  interface cap_series
    module procedure cap_series_vector
    module procedure cap_series_matrix
  end interface

  interface return_series
    module procedure return_series_vector
    module procedure return_series_matrix
  end interface

  interface return_convert
    module procedure return_convert_vector
    module procedure return_convert_matrix
  end interface

  interface capser
    module procedure cap_series_vector
    module procedure cap_series_matrix
  end interface

  interface returnseries
    module procedure return_series_vector
    module procedure return_series_matrix
  end interface

  interface returnconvert
    module procedure return_convert_vector
    module procedure return_convert_matrix
  end interface

  interface trdbinary
    module procedure trend_binary
    module procedure trend_binary_matrix
  end interface

  interface trdbilson
    module procedure trend_bilson
    module procedure trend_bilson_matrix
  end interface

  interface trdes
    module procedure trend_exponential
    module procedure trend_exponential_matrix
  end interface

  interface trdhp
    module procedure trend_hodrick_prescott
    module procedure trend_hodrick_prescott_matrix
  end interface

  interface trdsma
    module procedure trend_simple_moving_average
    module procedure trend_simple_moving_average_matrix
  end interface

  interface trdwma
    module procedure trend_weighted_moving_average
    module procedure trend_weighted_moving_average_matrix
  end interface

contains

  pure function cap_series_vector(y, minimum, maximum) result(out)
    real(dp), intent(in) :: y(:), minimum, maximum
    real(dp) :: out(size(y))
    out = min(max(y, minimum), maximum)
  end function cap_series_vector

  pure function cap_series_matrix(y, minimum, maximum) result(out)
    real(dp), intent(in) :: y(:, :), minimum, maximum
    real(dp) :: out(size(y, 1), size(y, 2))
    out = min(max(y, minimum), maximum)
  end function cap_series_matrix

  function return_series_vector(y, method, percentage, trim, compound, status) result(ret)
    real(dp), intent(in) :: y(:)
    integer, intent(in), optional :: method
    logical, intent(in), optional :: percentage, trim, compound
    integer, intent(out), optional :: status
    real(dp), allocatable :: ret(:)
    integer :: meth, i, start_index, n
    logical :: pct, do_trim, do_compound
    real(dp) :: scale, nan_value

    n = size(y)
    meth = returns_continuous
    if (present(method)) meth = method
    pct = .true.
    if (present(percentage)) pct = percentage
    do_trim = .false.
    if (present(trim)) do_trim = trim
    do_compound = .false.
    if (present(compound)) do_compound = compound
    start_index = merge(2, 1, do_trim)
    allocate(ret(max(0, n - start_index + 1)))
    if (n < 1 .or. (meth /= returns_continuous .and. meth /= returns_discrete)) then
      if (present(status)) status = frapo_invalid_input
      return
    end if
    if (meth == returns_continuous .and. any(y <= 0.0_dp)) then
      ret = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(status)) status = frapo_invalid_input
      return
    end if
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    block
      real(dp), allocatable :: full(:)
      allocate(full(n))
      full(1) = nan_value
      do i = 2, n
        if (meth == returns_continuous) then
          full(i) = log(y(i)) - log(y(i - 1))
        else
          full(i) = (y(i) - y(i - 1)) / y(i - 1)
        end if
      end do
      if (do_compound) then
        full(1) = 0.0_dp
        if (meth == returns_continuous) then
          do i = 2, n
            full(i) = full(i - 1) + full(i)
          end do
        else
          do i = 2, n
            full(i) = (1.0_dp + full(i - 1)) * (1.0_dp + full(i)) - 1.0_dp
          end do
        end if
      end if
      scale = merge(100.0_dp, 1.0_dp, pct)
      ret = scale * full(start_index:n)
    end block
    if (present(status)) status = frapo_ok
  end function return_series_vector

  function return_series_matrix(y, method, percentage, trim, compound, status) result(ret)
    real(dp), intent(in) :: y(:, :)
    integer, intent(in), optional :: method
    logical, intent(in), optional :: percentage, trim, compound
    integer, intent(out), optional :: status
    real(dp), allocatable :: ret(:, :), column(:)
    integer :: j, istat = frapo_ok, nout

    column = return_series_vector(y(:, 1), method, percentage, trim, compound, istat)
    nout = size(column)
    allocate(ret(nout, size(y, 2)))
    ret(:, 1) = column
    do j = 2, size(y, 2)
      column = return_series_vector(y(:, j), method, percentage, trim, compound, istat)
      ret(:, j) = column
    end do
    if (present(status)) status = istat
  end function return_series_matrix

  pure function return_convert_vector(y, direction, percentage) result(r)
    real(dp), intent(in) :: y(:)
    integer, intent(in), optional :: direction
    logical, intent(in), optional :: percentage
    real(dp) :: r(size(y))
    integer :: dir
    logical :: pct

    dir = convert_cont_to_disc
    if (present(direction)) dir = direction
    pct = .true.
    if (present(percentage)) pct = percentage
    r = y
    if (pct) r = r / 100.0_dp
    if (dir == convert_cont_to_disc) then
      r = exp(r) - 1.0_dp
    else
      r = log(1.0_dp + r)
    end if
    if (pct) r = 100.0_dp * r
  end function return_convert_vector

  pure function return_convert_matrix(y, direction, percentage) result(r)
    real(dp), intent(in) :: y(:, :)
    integer, intent(in), optional :: direction
    logical, intent(in), optional :: percentage
    real(dp) :: r(size(y, 1), size(y, 2))
    integer :: j
    do j = 1, size(y, 2)
      r(:, j) = return_convert_vector(y(:, j), direction, percentage)
    end do
  end function return_convert_matrix

  pure function trend_binary(y) result(trend)
    real(dp), intent(in) :: y(:)
    real(dp) :: trend(size(y))
    real(dp), parameter :: pi = acos(-1.0_dp)
    trend = sign(1.0_dp, y) * min(abs((4.0_dp / pi) * atan(y)), 1.0_dp)
    where (abs(y) <= tiny(1.0_dp)) trend = 0.0_dp
  end function trend_binary

  pure function trend_bilson(y, exponent) result(trend)
    real(dp), intent(in) :: y(:), exponent
    real(dp) :: trend(size(y))
    trend = sign(1.0_dp, y) * abs(y)**(1.0_dp - abs(y)**exponent)
    where (abs(y) <= tiny(1.0_dp)) trend = 0.0_dp
  end function trend_bilson

  function trend_exponential(y, lambda, initial, status) result(trend)
    real(dp), intent(in) :: y(:), lambda
    real(dp), intent(in), optional :: initial
    integer, intent(out), optional :: status
    real(dp), allocatable :: trend(:)
    real(dp) :: prior
    integer :: i

    allocate(trend(size(y)))
    if (lambda < 0.0_dp .or. lambda > 1.0_dp) then
      trend = 0.0_dp
      if (present(status)) status = frapo_invalid_input
      return
    end if
    prior = 0.0_dp
    if (present(initial)) prior = initial
    do i = 1, size(y)
      trend(i) = lambda * y(i) + (1.0_dp - lambda) * prior
      prior = trend(i)
    end do
    if (present(status)) status = frapo_ok
  end function trend_exponential

  function trend_hodrick_prescott(y, lambda, status) result(trend)
    real(dp), intent(in) :: y(:), lambda
    integer, intent(out), optional :: status
    real(dp), allocatable :: trend(:), system(:, :), solution(:)
    integer :: n, i, istat

    n = size(y)
    allocate(system(n, n))
    system = 0.0_dp
    do i = 1, n
      system(i, i) = 1.0_dp
    end do
    if (n >= 3) then
      do i = 1, n - 2
        system(i, i) = system(i, i) + lambda
        system(i, i + 1) = system(i, i + 1) - 2.0_dp * lambda
        system(i, i + 2) = system(i, i + 2) + lambda
        system(i + 1, i) = system(i + 1, i) - 2.0_dp * lambda
        system(i + 1, i + 1) = system(i + 1, i + 1) + 4.0_dp * lambda
        system(i + 1, i + 2) = system(i + 1, i + 2) - 2.0_dp * lambda
        system(i + 2, i) = system(i + 2, i) + lambda
        system(i + 2, i + 1) = system(i + 2, i + 1) - 2.0_dp * lambda
        system(i + 2, i + 2) = system(i + 2, i + 2) + lambda
      end do
    end if
    call solve_linear(system, y, solution, istat)
    call move_alloc(solution, trend)
    if (present(status)) status = istat
  end function trend_hodrick_prescott

  function trend_simple_moving_average(y, periods, trim, status) result(trend)
    real(dp), intent(in) :: y(:)
    integer, intent(in) :: periods
    logical, intent(in), optional :: trim
    integer, intent(out), optional :: status
    real(dp), allocatable :: trend(:), full(:)
    real(dp) :: nan_value
    integer :: k, i
    logical :: do_trim

    k = abs(periods)
    do_trim = .false.
    if (present(trim)) do_trim = trim
    if (k < 1 .or. k > size(y)) then
      allocate(trend(0))
      if (present(status)) status = frapo_invalid_input
      return
    end if
    allocate(full(size(y)))
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    full(1:k - 1) = nan_value
    do i = k, size(y)
      full(i) = sum(y(i - k + 1:i)) / real(k, dp)
    end do
    if (do_trim) then
      trend = full(k:)
    else
      trend = full
    end if
    if (present(status)) status = frapo_ok
  end function trend_simple_moving_average

  function trend_weighted_moving_average(y, weights, trim, status) result(trend)
    real(dp), intent(in) :: y(:), weights(:)
    logical, intent(in), optional :: trim
    integer, intent(out), optional :: status
    real(dp), allocatable :: trend(:), full(:)
    real(dp) :: nan_value
    integer :: k, i, j
    logical :: do_trim

    k = size(weights)
    do_trim = .false.
    if (present(trim)) do_trim = trim
    if (k < 1 .or. k > size(y)) then
      allocate(trend(0))
      if (present(status)) status = frapo_invalid_input
      return
    end if
    allocate(full(size(y)))
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    full(1:k - 1) = nan_value
    do i = k, size(y)
      full(i) = 0.0_dp
      do j = 1, k
        full(i) = full(i) + weights(j) * y(i - j + 1)
      end do
    end do
    if (do_trim) then
      trend = full(k:)
    else
      trend = full
    end if
    if (present(status)) status = frapo_ok
  end function trend_weighted_moving_average

  pure function trend_binary_matrix(y) result(trend)
    real(dp), intent(in) :: y(:, :)
    real(dp) :: trend(size(y, 1), size(y, 2))
    real(dp), parameter :: pi = acos(-1.0_dp)
    trend = sign(1.0_dp, y) * min(abs((4.0_dp / pi) * atan(y)), 1.0_dp)
    where (abs(y) <= tiny(1.0_dp)) trend = 0.0_dp
  end function trend_binary_matrix

  pure function trend_bilson_matrix(y, exponent) result(trend)
    real(dp), intent(in) :: y(:, :), exponent
    real(dp) :: trend(size(y, 1), size(y, 2))
    trend = sign(1.0_dp, y) * abs(y)**(1.0_dp - abs(y)**exponent)
    where (abs(y) <= tiny(1.0_dp)) trend = 0.0_dp
  end function trend_bilson_matrix

  function trend_exponential_matrix(y, lambda, initial, status) result(trend)
    real(dp), intent(in) :: y(:, :), lambda
    real(dp), intent(in), optional :: initial(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: trend(:, :), column(:)
    integer :: j, istat = frapo_ok

    allocate(trend(size(y, 1), size(y, 2)))
    if (present(initial)) then
      if (size(initial) /= size(y, 2)) then
        trend = 0.0_dp
        if (present(status)) status = frapo_invalid_input
        return
      end if
    end if
    do j = 1, size(y, 2)
      if (present(initial)) then
        column = trend_exponential(y(:, j), lambda, initial(j), istat)
      else
        column = trend_exponential(y(:, j), lambda, status=istat)
      end if
      trend(:, j) = column
    end do
    if (present(status)) status = istat
  end function trend_exponential_matrix

  function trend_hodrick_prescott_matrix(y, lambda, status) result(trend)
    real(dp), intent(in) :: y(:, :), lambda
    integer, intent(out), optional :: status
    real(dp), allocatable :: trend(:, :), column(:)
    integer :: j, istat = frapo_ok

    allocate(trend(size(y, 1), size(y, 2)))
    do j = 1, size(y, 2)
      column = trend_hodrick_prescott(y(:, j), lambda, istat)
      trend(:, j) = column
    end do
    if (present(status)) status = istat
  end function trend_hodrick_prescott_matrix

  function trend_simple_moving_average_matrix(y, periods, trim, status) result(trend)
    real(dp), intent(in) :: y(:, :)
    integer, intent(in) :: periods
    logical, intent(in), optional :: trim
    integer, intent(out), optional :: status
    real(dp), allocatable :: trend(:, :), column(:)
    integer :: j, istat = frapo_ok, nout

    column = trend_simple_moving_average(y(:, 1), periods, trim, istat)
    nout = size(column)
    allocate(trend(nout, size(y, 2)))
    trend(:, 1) = column
    do j = 2, size(y, 2)
      column = trend_simple_moving_average(y(:, j), periods, trim, istat)
      trend(:, j) = column
    end do
    if (present(status)) status = istat
  end function trend_simple_moving_average_matrix

  function trend_weighted_moving_average_matrix(y, weights, trim, status) result(trend)
    real(dp), intent(in) :: y(:, :), weights(:)
    logical, intent(in), optional :: trim
    integer, intent(out), optional :: status
    real(dp), allocatable :: trend(:, :), column(:)
    integer :: j, istat = frapo_ok, nout

    column = trend_weighted_moving_average(y(:, 1), weights, trim, istat)
    nout = size(column)
    allocate(trend(nout, size(y, 2)))
    trend(:, 1) = column
    do j = 2, size(y, 2)
      column = trend_weighted_moving_average(y(:, j), weights, trim, istat)
      trend(:, j) = column
    end do
    if (present(status)) status = istat
  end function trend_weighted_moving_average_matrix
end module frapo_series
