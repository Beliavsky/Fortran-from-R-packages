! SPDX-License-Identifier: GPL-3.0-only
module mass_nonlinear
  use rrcov_kinds, only : dp
  use mass_types, only : mass_success, mass_invalid_argument, mass_dimension_error
  use mass_math, only : least_squares
  implicit none
  private
  public :: negative_exponential, negexp_initial
contains

  pure elemental function negative_exponential(x, b0, b1, theta) result(value)
    real(dp), intent(in) :: x, b0, b1, theta
    real(dp) :: value
    value = b0 + b1 * exp(-x / theta)
  end function negative_exponential

  subroutine negexp_initial(x, y, parameters, status)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(out) :: parameters(3)
    integer, intent(out) :: status
    real(dp), allocatable :: design(:, :), beta(:), residuals(:)
    real(dp), allocatable :: x0(:), values(:), design2(:, :)
    real(dp) :: mean_x, xh, lo, hi, theta
    integer :: rank, st
    if (size(x) /= size(y) .or. size(x) < 3) then
      parameters = 0.0_dp
      status = mass_dimension_error
      return
    end if
    if (maxval(x) <= minval(x)) then
      parameters = 0.0_dp
      status = mass_invalid_argument
      return
    end if
    mean_x = sum(x) / real(size(x), dp)
    allocate(design(size(x), 3))
    design(:, 1) = 1.0_dp
    design(:, 2) = x - mean_x
    design(:, 3) = -(x - mean_x)**2 / 2.0_dp
    call least_squares(design, y, beta, residuals, rank, st)
    if (st /= mass_success .or. abs(beta(3)) <= sqrt(epsilon(1.0_dp))) then
      parameters = 0.0_dp
      status = mass_invalid_argument
      return
    end if
    lo = minval(x)
    hi = maxval(x)
    xh = mean_x + beta(2) / beta(3)
    if ((xh - lo) * (xh - hi) < 0.0_dp) then
      if (xh - lo > hi - xh) then
        hi = xh
      else
        lo = xh
      end if
    end if
    allocate(x0(3), values(3))
    x0 = [lo, 0.5_dp * (lo + hi), hi]
    values = beta(1) + beta(2) * (x0 - mean_x) - &
      beta(3) * (x0 - mean_x)**2 / 2.0_dp
    if ((values(2) - values(1)) * (values(3) - values(2)) <= 0.0_dp) then
      parameters = 0.0_dp
      status = mass_invalid_argument
      return
    end if
    theta = (x0(2) - x0(1)) / log((values(2) - values(1)) / &
      (values(3) - values(2)))
    if (.not. (theta > 0.0_dp)) then
      parameters = 0.0_dp
      status = mass_invalid_argument
      return
    end if
    allocate(design2(size(x), 2))
    design2(:, 1) = 1.0_dp
    design2(:, 2) = exp(-x / theta)
    call least_squares(design2, y, beta, residuals, rank, st)
    parameters = [beta(1), beta(2), theta]
    status = st
  end subroutine negexp_initial

end module mass_nonlinear
