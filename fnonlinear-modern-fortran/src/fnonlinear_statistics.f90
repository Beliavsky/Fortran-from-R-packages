! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
module fnonlinear_statistics
  use chaos_kinds, only : dp
  use chaos_embedding, only : delay_embed
  use chaos_metrics, only : average_mutual_information, correlation_integral, &
    correlation_dimension_curve, recurrence_distance_matrix, space_time_separation
  use chaos_neighbors, only : false_nearest_curve, find_k_nearests, &
    lyapunov_stretching, lyapunov_linear_fit
  implicit none
  private
  public :: mutual_information_curve, false_nearest_neighbors
  public :: recurrence_matrix, recurrence_distance_matrix
  public :: space_time_separation, lyapunov_stretching, lyapunov_linear_fit
  public :: correlation_integral, correlation_dimension_curve, find_k_nearests
contains
  subroutine mutual_information_curve(x, partitions, lag_max, values, status)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: partitions, lag_max
    real(dp), allocatable, intent(out) :: values(:)
    integer, intent(out) :: status
    call average_mutual_information(x, partitions, lag_max, values, status)
  end subroutine mutual_information_curve

  subroutine false_nearest_neighbors(x, max_dimension, delay, theiler, escape_factor, radius, &
      fractions, totals, status, search_method)
    real(dp), intent(in) :: x(:), escape_factor, radius
    integer, intent(in) :: max_dimension, delay, theiler
    real(dp), allocatable, intent(out) :: fractions(:)
    integer, allocatable, intent(out) :: totals(:)
    integer, intent(out) :: status
    character(len=*), intent(in), optional :: search_method
    call false_nearest_curve(x, max_dimension, delay, theiler, escape_factor, radius, &
      fractions, totals, status, search_method=search_method)
  end subroutine false_nearest_neighbors

  subroutine recurrence_matrix(x, dimension, delay, end_time, epsilon, recurrence, status, step)
    real(dp), intent(in) :: x(:), epsilon
    integer, intent(in) :: dimension, delay, end_time
    logical, allocatable, intent(out) :: recurrence(:, :)
    integer, intent(out) :: status
    integer, intent(in), optional :: step
    real(dp), allocatable :: embedded(:, :)
    integer :: embed_status, npoint, stride, i, j, ii, jj
    real(dp) :: epsilon2

    if (epsilon <= 0.0_dp .or. end_time < 1) then
      allocate(recurrence(0, 0))
      status = 1
      return
    end if
    call delay_embed(x, dimension, delay, embedded, embed_status)
    if (embed_status /= 0) then
      allocate(recurrence(0, 0))
      status = 2
      return
    end if
    stride = 1
    if (present(step)) stride = max(1, step)
    npoint = min(end_time, size(embedded, 1))
    npoint = (npoint - 1) / stride + 1
    allocate(recurrence(npoint, npoint))
    epsilon2 = epsilon**2
    do i = 1, npoint
      ii = 1 + (i - 1) * stride
      do j = 1, npoint
        jj = 1 + (j - 1) * stride
        recurrence(i, j) = sum((embedded(ii, :) - embedded(jj, :))**2) < epsilon2
      end do
    end do
    status = 0
  end subroutine recurrence_matrix
end module fnonlinear_statistics
