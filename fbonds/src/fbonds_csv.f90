! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! Based on fBonds, copyright its original authors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License, version 2 or later.
module fbonds_csv
   use fbonds_kinds, only : dp
   implicit none
   private
   public :: read_maturity_rate_csv
contains
   subroutine read_maturity_rate_csv(filename, maturity, rate, status)
      character(len=*), intent(in) :: filename
      real(dp), allocatable, intent(out) :: maturity(:), rate(:)
      integer, intent(out) :: status
      character(len=2048) :: line
      real(dp) :: x, y
      integer :: unit, ios, count, i

      status = 0
      count = 0
      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         allocate(maturity(0), rate(0))
         status = ios
         return
      end if
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         call commas_to_spaces(line)
         read(line, *, iostat=ios) x, y
         if (ios == 0) count = count + 1
      end do
      rewind(unit)
      allocate(maturity(count), rate(count))
      i = 0
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         call commas_to_spaces(line)
         read(line, *, iostat=ios) x, y
         if (ios == 0) then
            i = i + 1
            maturity(i) = x
            rate(i) = y
         end if
      end do
      close(unit)
      if (count == 0) status = -1
   end subroutine read_maturity_rate_csv

   pure subroutine commas_to_spaces(line)
      character(len=*), intent(inout) :: line
      integer :: i
      do i = 1, len_trim(line)
         if (line(i:i) == ',') line(i:i) = ' '
      end do
   end subroutine commas_to_spaces
end module fbonds_csv
