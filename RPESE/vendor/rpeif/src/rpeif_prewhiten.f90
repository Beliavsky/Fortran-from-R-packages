! SPDX-License-Identifier: GPL-3.0-or-later
module rpeif_prewhiten
  use rpeif_kinds, only : dp
  implicit none
  private
  public :: ar_prewhiten
contains
  subroutine ar_prewhiten(x, order, residuals, coefficients, intercept, status)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: order
    real(dp), allocatable, intent(out) :: residuals(:)
    real(dp), allocatable, intent(out), optional :: coefficients(:)
    real(dp), intent(out), optional :: intercept
    integer, intent(out), optional :: status
    real(dp), allocatable :: normal_matrix(:, :), rhs(:), beta(:), row(:)
    integer :: n, p, t, j, solve_status

    n = size(x)
    p = order
    allocate(residuals(n))
    residuals = 0.0_dp
    if (p < 1 .or. n <= p + 1) then
      residuals = x
      if (present(coefficients)) allocate(coefficients(0))
      if (present(intercept)) intercept = 0.0_dp
      if (present(status)) status = 1
      return
    end if

    allocate(normal_matrix(p + 1, p + 1), rhs(p + 1), beta(p + 1), row(p + 1))
    normal_matrix = 0.0_dp
    rhs = 0.0_dp
    do t = p + 1, n
      row(1) = 1.0_dp
      do j = 1, p
        row(j + 1) = x(t - j)
      end do
      normal_matrix = normal_matrix + spread(row, 2, p + 1) * spread(row, 1, p + 1)
      rhs = rhs + row * x(t)
    end do
    call solve_linear_system(normal_matrix, rhs, beta, solve_status)
    if (solve_status /= 0) then
      residuals = x - sum(x) / real(n, dp)
      if (present(coefficients)) then
        allocate(coefficients(p))
        coefficients = 0.0_dp
      end if
      if (present(intercept)) intercept = sum(x) / real(n, dp)
      if (present(status)) status = 2
      return
    end if

    residuals(1:p) = 0.0_dp
    do t = p + 1, n
      residuals(t) = x(t) - beta(1)
      do j = 1, p
        residuals(t) = residuals(t) - beta(j + 1) * x(t - j)
      end do
    end do
    if (present(coefficients)) then
      allocate(coefficients(p))
      coefficients = beta(2:)
    end if
    if (present(intercept)) intercept = beta(1)
    if (present(status)) status = 0
  end subroutine ar_prewhiten

  subroutine solve_linear_system(a, b, x, status)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), allocatable :: aug(:, :), tmp_row(:)
    real(dp) :: factor, pivot_value
    integer :: n, i, j, pivot

    n = size(b)
    status = 0
    allocate(aug(n, n + 1), tmp_row(n + 1))
    aug(:, 1:n) = a
    aug(:, n + 1) = b

    do i = 1, n
      pivot = i - 1 + maxloc(abs(aug(i:n, i)), dim=1)
      pivot_value = abs(aug(pivot, i))
      if (pivot_value <= 1.0e-14_dp) then
        status = 1
        x = 0.0_dp
        return
      end if
      if (pivot /= i) then
        tmp_row = aug(i, :)
        aug(i, :) = aug(pivot, :)
        aug(pivot, :) = tmp_row
      end if
      aug(i, :) = aug(i, :) / aug(i, i)
      do j = 1, n
        if (j == i) cycle
        factor = aug(j, i)
        aug(j, :) = aug(j, :) - factor * aug(i, :)
      end do
    end do
    x = aug(:, n + 1)
  end subroutine solve_linear_system
end module rpeif_prewhiten
