! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of tseriesChaos and is distributed
! under the GNU General Public License version 2 only.
module chaos_embedding
  use chaos_kinds, only : dp
  implicit none
  private
  public :: delay_embed, delay_embed_lags, delay_embed_matrix
contains
  subroutine delay_embed(x, m, d, embedded, status)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: m, d
    real(dp), allocatable, intent(out) :: embedded(:, :)
    integer, intent(out) :: status
    integer, allocatable :: lags(:)
    integer :: j

    if (m < 1 .or. d < 1) then
      allocate(embedded(0, 0))
      status = 1
      return
    end if
    allocate(lags(m))
    do j = 1, m
      lags(j) = (j - 1) * d
    end do
    call delay_embed_lags(x, lags, embedded, status)
  end subroutine delay_embed

  subroutine delay_embed_lags(x, lags, embedded, status)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lags(:)
    real(dp), allocatable, intent(out) :: embedded(:, :)
    integer, intent(out) :: status
    integer :: nrow, i, j, maxlag

    if (size(lags) == 0 .or. any(lags < 0)) then
      allocate(embedded(0, 0))
      status = 1
      return
    end if
    maxlag = maxval(lags)
    nrow = size(x) - maxlag
    if (nrow <= 0) then
      allocate(embedded(0, size(lags)))
      status = 2
      return
    end if
    allocate(embedded(nrow, size(lags)))
    do j = 1, size(lags)
      do i = 1, nrow
        embedded(i, j) = x(i + lags(j))
      end do
    end do
    status = 0
  end subroutine delay_embed_lags

  subroutine delay_embed_matrix(x, lags, embedded, status)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: lags(:)
    real(dp), allocatable, intent(out) :: embedded(:, :)
    integer, intent(out) :: status
    integer :: nrow, maxlag, i, j, k, col

    if (size(lags) == 0 .or. any(lags < 0)) then
      allocate(embedded(0, 0))
      status = 1
      return
    end if
    maxlag = maxval(lags)
    nrow = size(x, 1) - maxlag
    if (nrow <= 0) then
      allocate(embedded(0, size(x, 2) * size(lags)))
      status = 2
      return
    end if
    allocate(embedded(nrow, size(x, 2) * size(lags)))
    col = 0
    do j = 1, size(lags)
      do k = 1, size(x, 2)
        col = col + 1
        do i = 1, nrow
          embedded(i, col) = x(i + lags(j), k)
        end do
      end do
    end do
    status = 0
  end subroutine delay_embed_matrix
end module chaos_embedding
