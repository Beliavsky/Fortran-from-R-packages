! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
module garchx_csv
   use garchx_kinds, only : dp
   implicit none
   private
   public :: read_dated_numeric_csv, parse_integer_list
contains
   subroutine replace_commas(line)
      character(len=*), intent(inout) :: line
      integer :: i
      do i = 1, len_trim(line)
         if (line(i:i) == ',') line(i:i) = ' '
      end do
   end subroutine replace_commas

   integer function comma_count(line) result(count_value)
      character(len=*), intent(in) :: line
      integer :: i
      count_value = 0
      do i = 1, len_trim(line)
         if (line(i:i) == ',') count_value = count_value+1
      end do
   end function comma_count

   subroutine read_dated_numeric_csv(filename, dates, values, status)
      character(len=*), intent(in) :: filename
      character(len=64), allocatable, intent(out) :: dates(:)
      real(dp), allocatable, intent(out) :: values(:, :)
      integer, intent(out) :: status
      integer :: unit, ios, nrow, ncol, i
      character(len=8192) :: line

      status = 0
      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         status = 1
         allocate(dates(0), values(0, 0))
         return
      end if
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) then
         status = 2
         close(unit)
         allocate(dates(0), values(0, 0))
         return
      end if
      ncol = comma_count(line)
      if (ncol < 1) then
         status = 3
         close(unit)
         allocate(dates(0), values(0, 0))
         return
      end if
      nrow = 0
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) > 0) nrow = nrow+1
      end do
      if (nrow < 1) then
         status = 4
         close(unit)
         allocate(dates(0), values(0, 0))
         return
      end if
      rewind(unit)
      read(unit, '(a)') line
      allocate(dates(nrow), values(nrow, ncol))
      i = 0
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         i = i+1
         call replace_commas(line)
         read(line, *, iostat=ios) dates(i), values(i, :)
         if (ios /= 0) then
            status = 5
            close(unit)
            return
         end if
      end do
      close(unit)
   end subroutine read_dated_numeric_csv

   subroutine parse_integer_list(text, values, status)
      character(len=*), intent(in) :: text
      integer, allocatable, intent(out) :: values(:)
      integer, intent(out) :: status
      character(len=:), allocatable :: work
      integer :: n, ios

      work = trim(adjustl(text))
      if (len_trim(work) == 0 .or. work == '-' .or. work == 'none' .or. work == 'NULL') then
         allocate(values(0))
         status = 0
         return
      end if
      n = comma_count(work)+1
      allocate(values(n))
      call replace_commas(work)
      read(work, *, iostat=ios) values
      if (ios /= 0 .or. any(values < 1)) then
         status = 1
         values = 0
      else
         status = 0
      end if
   end subroutine parse_integer_list
end module garchx_csv
