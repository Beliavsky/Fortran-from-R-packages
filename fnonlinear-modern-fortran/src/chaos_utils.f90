! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of fNonlinear/tseriesChaos computational routines and is distributed
! under the GNU General Public License version 2 or later.
module chaos_utils
  use chaos_kinds, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: embedding_length, scale_unit_interval, squared_embedding_distance
  public :: linear_regression, sort_pairs, quiet_nan
contains
  pure integer function embedding_length(n, m, d) result(nembed)
    integer, intent(in) :: n, m, d
    if (n < 0 .or. m < 1 .or. d < 1) then
      nembed = 0
    else
      nembed = n - (m - 1) * d
    end if
  end function embedding_length

  subroutine scale_unit_interval(x, y, span, status)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: y(:)
    real(dp), intent(out) :: span
    integer, intent(out) :: status
    real(dp) :: xmin, xmax

    allocate(y(size(x)))
    if (size(x) == 0) then
      span = 0.0_dp
      status = 1
      return
    end if
    xmin = minval(x)
    xmax = maxval(x)
    span = xmax - xmin
    if (span <= 0.0_dp) then
      y = 0.0_dp
      status = 2
    else
      y = (x - xmin) / span
      status = 0
    end if
  end subroutine scale_unit_interval

  pure real(dp) function squared_embedding_distance(x, ia, ib, m, d) result(dist2)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: ia, ib, m, d
    integer :: j
    dist2 = 0.0_dp
    do j = 0, m - 1
      dist2 = dist2 + (x(ia + j * d) - x(ib + j * d))**2
    end do
  end function squared_embedding_distance

  subroutine linear_regression(x, y, intercept, slope, status)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(out) :: intercept, slope
    integer, intent(out) :: status
    real(dp) :: xbar, ybar, sxx, sxy

    if (size(x) /= size(y) .or. size(x) < 2) then
      intercept = quiet_nan()
      slope = quiet_nan()
      status = 1
      return
    end if
    xbar = sum(x) / real(size(x), dp)
    ybar = sum(y) / real(size(y), dp)
    sxx = sum((x - xbar)**2)
    if (sxx <= tiny(1.0_dp)) then
      intercept = quiet_nan()
      slope = quiet_nan()
      status = 2
      return
    end if
    sxy = sum((x - xbar) * (y - ybar))
    slope = sxy / sxx
    intercept = ybar - slope * xbar
    status = 0
  end subroutine linear_regression

  subroutine sort_pairs(values, indices, n)
    real(dp), intent(inout) :: values(:)
    integer, intent(inout) :: indices(:)
    integer, intent(in) :: n
    integer :: i, j, itmp
    real(dp) :: vtmp

    do i = 2, n
      vtmp = values(i)
      itmp = indices(i)
      j = i - 1
      do while (j >= 1)
        if (values(j) < vtmp) exit
        if (.not. (values(j) > vtmp)) then
          if (indices(j) <= itmp) exit
        end if
        values(j + 1) = values(j)
        indices(j + 1) = indices(j)
        j = j - 1
      end do
      values(j + 1) = vtmp
      indices(j + 1) = itmp
    end do
  end subroutine sort_pairs

  pure real(dp) function quiet_nan() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function quiet_nan
end module chaos_utils
