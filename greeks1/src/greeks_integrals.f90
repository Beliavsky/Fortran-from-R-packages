! SPDX-License-Identifier: MIT
module greeks_integrals
  use greeks_kinds, only: dp
  implicit none
  private
  public :: row_cumsums, make_bm, calc_i, calc_i_1, calc_i_2, calc_i_3
  public :: calc_x, calc_log_x, calc_xw, calc_txw

contains

  subroutine row_cumsums(matrix)
    real(dp), intent(inout) :: matrix(:, :)
    integer :: j
    do j = 2, size(matrix, 2)
      matrix(:, j) = matrix(:, j - 1) + matrix(:, j)
    end do
  end subroutine row_cumsums

  subroutine make_bm(increments, paths, steps, w)
    real(dp), intent(in) :: increments(:)
    integer, intent(in) :: paths, steps
    real(dp), allocatable, intent(out) :: w(:, :)
    integer :: i, j, index_value

    allocate(w(paths, steps + 1))
    w(:, 1) = 0.0_dp
    index_value = 0
    do j = 2, steps + 1
      do i = 1, paths
        index_value = index_value + 1
        w(i, j) = w(i, j - 1) + increments(index_value)
      end do
    end do
  end subroutine make_bm

  function calc_i(x, dt) result(values)
    real(dp), intent(in) :: x(:, :), dt
    real(dp), allocatable :: values(:)
    integer :: ncols
    ncols = size(x, 2)
    allocate(values(size(x, 1)))
    values = 0.5_dp*(x(:, 1) + x(:, ncols))
    if (ncols > 2) values = values + sum(x(:, 2:ncols - 1), dim=2)
    values = values*dt
  end function calc_i

  function calc_i_1(x, dt) result(values)
    real(dp), intent(in) :: x(:, :), dt
    real(dp), allocatable :: values(:)
    integer :: j, steps
    steps = size(x, 2) - 1
    allocate(values(size(x, 1)))
    values = 0.5_dp*x(:, steps + 1)*real(steps, dp)*dt
    do j = 2, steps
      values = values + x(:, j)*real(j - 1, dp)*dt
    end do
    values = values*dt
  end function calc_i_1

  function calc_i_2(x, dt) result(values)
    real(dp), intent(in) :: x(:, :), dt
    real(dp), allocatable :: values(:)
    integer :: j, steps
    real(dp) :: time_value
    steps = size(x, 2) - 1
    allocate(values(size(x, 1)))
    time_value = real(steps, dp)*dt
    values = 0.5_dp*x(:, steps + 1)*time_value**2
    do j = 2, steps
      time_value = real(j - 1, dp)*dt
      values = values + x(:, j)*time_value**2
    end do
    values = values*dt
  end function calc_i_2

  function calc_i_3(x, dt) result(values)
    real(dp), intent(in) :: x(:, :), dt
    real(dp), allocatable :: values(:)
    integer :: j, steps
    real(dp) :: time_value
    steps = size(x, 2) - 1
    allocate(values(size(x, 1)))
    time_value = real(steps, dp)*dt
    values = 0.5_dp*x(:, steps + 1)*time_value**3
    do j = 2, steps
      time_value = real(j - 1, dp)*dt
      values = values + x(:, j)*time_value**3
    end do
    values = values*dt
  end function calc_i_3

  function calc_x(w, dt, volatility, drift) result(x)
    real(dp), intent(in) :: w(:, :), dt, volatility, drift
    real(dp), allocatable :: x(:, :)
    integer :: j
    allocate(x(size(w, 1), size(w, 2)))
    do j = 1, size(w, 2)
      x(:, j) = exp((drift - 0.5_dp*volatility**2)*real(j - 1, dp)*dt + &
        volatility*w(:, j))
    end do
  end function calc_x

  function calc_log_x(w, dt, volatility, drift) result(x)
    real(dp), intent(in) :: w(:, :), dt, volatility, drift
    real(dp), allocatable :: x(:, :)
    integer :: j
    allocate(x(size(w, 1), size(w, 2)))
    do j = 1, size(w, 2)
      x(:, j) = (drift - 0.5_dp*volatility**2)*real(j - 1, dp)*dt + &
        volatility*w(:, j)
    end do
  end function calc_log_x

  function calc_xw(x, w, dt) result(values)
    real(dp), intent(in) :: x(:, :), w(:, :), dt
    real(dp), allocatable :: values(:)
    integer :: ncols
    ncols = size(x, 2)
    allocate(values(size(x, 1)))
    values = 0.5_dp*x(:, ncols)*w(:, ncols)
    if (ncols > 2) values = values + sum(x(:, 2:ncols - 1)*w(:, 2:ncols - 1), dim=2)
    values = values*dt
  end function calc_xw

  function calc_txw(x, w, dt) result(values)
    real(dp), intent(in) :: x(:, :), w(:, :), dt
    real(dp), allocatable :: values(:)
    integer :: j, steps
    steps = size(x, 2) - 1
    allocate(values(size(x, 1)))
    values = 0.5_dp*x(:, steps + 1)*w(:, steps + 1)*real(steps, dp)*dt*dt
    do j = 2, steps
      values = values + x(:, j)*w(:, j)*real(j - 1, dp)*dt*dt
    end do
  end function calc_txw

end module greeks_integrals
