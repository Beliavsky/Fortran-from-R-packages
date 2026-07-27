! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
module garchx_utilities
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use garchx_kinds, only : dp
   implicit none
   private
   public :: glag, gdiff
   interface glag
      module procedure glag_vector, glag_matrix
   end interface glag
   interface gdiff
      module procedure gdiff_vector, gdiff_matrix
   end interface gdiff
contains
   subroutine glag_vector(x, k, result, status, pad, pad_value)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: result(:)
      integer, intent(out) :: status
      logical, intent(in), optional :: pad
      real(dp), intent(in), optional :: pad_value
      logical :: do_pad
      real(dp) :: fill
      if (k < 1 .or. k >= size(x)) then
         status = 1
         allocate(result(0))
         return
      end if
      do_pad = .true.
      if (present(pad)) do_pad = pad
      fill = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(pad_value)) fill = pad_value
      if (do_pad) then
         allocate(result(size(x)))
         result(1:k) = fill
         result(k+1:) = x(:size(x)-k)
      else
         allocate(result(size(x)-k))
         result = x(:size(x)-k)
      end if
      status = 0
   end subroutine glag_vector

   subroutine glag_matrix(x, k, result, status, pad, pad_value)
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: result(:, :)
      integer, intent(out) :: status
      logical, intent(in), optional :: pad
      real(dp), intent(in), optional :: pad_value
      logical :: do_pad
      real(dp) :: fill
      if (k < 1 .or. k >= size(x, 1)) then
         status = 1
         allocate(result(0, 0))
         return
      end if
      do_pad = .true.
      if (present(pad)) do_pad = pad
      fill = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(pad_value)) fill = pad_value
      if (do_pad) then
         allocate(result(size(x, 1), size(x, 2)))
         result(1:k, :) = fill
         result(k+1:, :) = x(:size(x, 1)-k, :)
      else
         allocate(result(size(x, 1)-k, size(x, 2)))
         result = x(:size(x, 1)-k, :)
      end if
      status = 0
   end subroutine glag_matrix

   subroutine gdiff_vector(x, lag, result, status, pad, pad_value)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: lag
      real(dp), allocatable, intent(out) :: result(:)
      integer, intent(out) :: status
      logical, intent(in), optional :: pad
      real(dp), intent(in), optional :: pad_value
      real(dp), allocatable :: lagged(:)
      logical :: do_pad
      if (lag < 1 .or. lag >= size(x)) then
         status = 1
         allocate(result(0))
         return
      end if
      do_pad = .true.
      if (present(pad)) do_pad = pad
      call glag_vector(x, lag, lagged, status, do_pad, pad_value)
      if (status /= 0) then
         allocate(result(0))
         return
      end if
      if (do_pad) then
         allocate(result(size(x)))
         result = x-lagged
         if (present(pad_value)) result(1:lag) = pad_value
      else
         allocate(result(size(x)-lag))
         result = x(lag+1:)-lagged
      end if
   end subroutine gdiff_vector

   subroutine gdiff_matrix(x, lag, result, status, pad, pad_value)
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: lag
      real(dp), allocatable, intent(out) :: result(:, :)
      integer, intent(out) :: status
      logical, intent(in), optional :: pad
      real(dp), intent(in), optional :: pad_value
      real(dp), allocatable :: lagged(:, :)
      logical :: do_pad
      if (lag < 1 .or. lag >= size(x, 1)) then
         status = 1
         allocate(result(0, 0))
         return
      end if
      do_pad = .true.
      if (present(pad)) do_pad = pad
      call glag_matrix(x, lag, lagged, status, do_pad, pad_value)
      if (status /= 0) then
         allocate(result(0, 0))
         return
      end if
      if (do_pad) then
         allocate(result(size(x, 1), size(x, 2)))
         result = x-lagged
         if (present(pad_value)) result(1:lag, :) = pad_value
      else
         allocate(result(size(x, 1)-lag, size(x, 2)))
         result = x(lag+1:, :)-lagged
      end if
   end subroutine gdiff_matrix
end module garchx_utilities
