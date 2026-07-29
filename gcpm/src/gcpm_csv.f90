! SPDX-License-Identifier: GPL-2.0-only
!
! Modern Fortran translation of computational methods from GCPM 1.2.2.
! Original software copyright (C) 2015 Kevin Jakob and Dr. Matthias Fischer.
! Fortran translation copyright (C) 2026.
module gcpm_csv
   use gcpm_kinds, only: dp
   use gcpm_types, only: credit_portfolio, default_bernoulli, default_poisson
   use gcpm_portfolio, only: allocate_portfolio
   implicit none
   private

   integer, parameter :: max_fields = 256
   integer, parameter :: field_len = 256
   integer, parameter :: line_len = 65536

   public :: read_gcpm_portfolio

contains

   subroutine read_gcpm_portfolio(filename, portfolio, status, message, delimiter, decimal_comma)
      character(len=*), intent(in) :: filename
      type(credit_portfolio), intent(out) :: portfolio
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      character(len=1), intent(in), optional :: delimiter
      logical, intent(in), optional :: decimal_comma

      character(len=1) :: sep
      logical :: comma_decimal
      character(len=line_len) :: line
      character(len=field_len) :: fields(max_fields)
      integer :: unit, ios, n_fields, n_rows, row, k, n_sectors

      status = 0
      message = ''
      sep = ';'
      if (present(delimiter)) sep = delimiter
      comma_decimal = sep == ';'
      if (present(decimal_comma)) comma_decimal = decimal_comma

      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         status = 1
         message = 'cannot open portfolio file: ' // trim(filename)
         return
      end if

      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) then
         close(unit)
         status = 2
         message = 'cannot read portfolio header'
         return
      end if
      call split_delimited(trim(line), sep, fields, n_fields)
      if (n_fields < 9) then
         close(unit)
         status = 3
         message = 'portfolio file must contain the eight GCPM columns plus sector weights'
         return
      end if
      n_sectors = n_fields - 8

      n_rows = 0
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) > 0) n_rows = n_rows + 1
      end do
      if (n_rows == 0) then
         close(unit)
         status = 4
         message = 'portfolio file contains no data rows'
         return
      end if

      rewind(unit)
      read(unit, '(a)') line
      call allocate_portfolio(portfolio, n_rows, n_sectors)
      call split_delimited(trim(line), sep, fields, n_fields)
      do k = 1, n_sectors
         portfolio%sector_name(k) = trim(fields(k + 8))
      end do

      row = 0
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         row = row + 1
         call split_delimited(trim(line), sep, fields, n_fields)
         if (n_fields /= n_sectors + 8) then
            close(unit)
            status = 5
            message = 'inconsistent field count at portfolio data row'
            return
         end if

         call read_integer_field(fields(1), portfolio%number(row), ios)
         if (ios /= 0) then
            status = 6
            message = 'invalid Number field at data row'
            close(unit)
            return
         end if
         portfolio%name(row) = trim(fields(2))
         call read_real_field(fields(5), comma_decimal, portfolio%ead(row), ios)
         if (ios /= 0) then
            status = 7
            message = 'invalid EAD field at data row'
            close(unit)
            return
         end if
         call read_real_field(fields(6), comma_decimal, portfolio%lgd(row), ios)
         if (ios /= 0) then
            status = 8
            message = 'invalid LGD field at data row'
            close(unit)
            return
         end if
         call read_real_field(fields(7), comma_decimal, portfolio%pd(row), ios)
         if (ios /= 0) then
            status = 9
            message = 'invalid PD field at data row'
            close(unit)
            return
         end if

         select case (trim(adjustl(fields(8))))
         case ('Bernoulli', 'bernoulli', 'BERNOULLI')
            portfolio%default_kind(row) = default_bernoulli
         case ('Poisson', 'poisson', 'POISSON')
            portfolio%default_kind(row) = default_poisson
         case default
            status = 10
            message = 'Default must be Bernoulli or Poisson'
            close(unit)
            return
         end select

         do k = 1, n_sectors
            call read_real_field(fields(k + 8), comma_decimal, portfolio%weight(row, k), ios)
            if (ios /= 0) then
               status = 11
               message = 'invalid sector-weight field at data row'
               close(unit)
               return
            end if
         end do
      end do
      close(unit)
   end subroutine read_gcpm_portfolio

   subroutine split_delimited(line, delimiter, fields, n_fields)
      character(len=*), intent(in) :: line
      character(len=1), intent(in) :: delimiter
      character(len=field_len), intent(out) :: fields(max_fields)
      integer, intent(out) :: n_fields
      integer :: i, start, length_line

      fields = ''
      n_fields = 0
      start = 1
      length_line = len_trim(line)
      do i = 1, length_line
         if (line(i:i) == delimiter) then
            n_fields = n_fields + 1
            if (n_fields <= max_fields) then
               if (i > start) fields(n_fields) = adjustl(line(start:i - 1))
            end if
            start = i + 1
         end if
      end do
      n_fields = n_fields + 1
      if (n_fields <= max_fields .and. start <= length_line) then
         fields(n_fields) = adjustl(line(start:length_line))
      end if
   end subroutine split_delimited

   subroutine read_real_field(field, decimal_comma, value, ios)
      character(len=*), intent(in) :: field
      logical, intent(in) :: decimal_comma
      real(dp), intent(out) :: value
      integer, intent(out) :: ios
      character(len=field_len) :: work
      integer :: i

      work = adjustl(field)
      if (decimal_comma) then
         do i = 1, len_trim(work)
            if (work(i:i) == ',') work(i:i) = '.'
         end do
      end if
      read(work, *, iostat=ios) value
   end subroutine read_real_field

   subroutine read_integer_field(field, value, ios)
      character(len=*), intent(in) :: field
      integer, intent(out) :: value
      integer, intent(out) :: ios

      read(field, *, iostat=ios) value
   end subroutine read_integer_field

end module gcpm_csv
