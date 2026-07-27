! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
program fit_csv
   use rmgarch
   implicit none

   character(len=1024) :: filename
   real(dp), allocatable :: data(:,:)
   type(multivariate_garch_fit) :: fit
   integer :: i

   if (command_argument_count() < 1) then
      print '(a)', 'usage: fpm run --example fit_csv -- FILE.csv'
      error stop 1
   end if
   call get_command_argument(1,filename)
   call read_csv_with_header(trim(filename),data)
   fit = fit_two_step_dcc(data,max_iterations=500)

   print '(a,i0)', 'observations: ', size(data,1)
   print '(a,i0)', 'assets: ', size(data,2)
   do i = 1, size(data,2)
      print '(a,i0,a,f10.6,a,f10.6,a,f10.6)', 'margin ',i,': alpha=',fit%margins(i)%alpha, &
         ' beta=',fit%margins(i)%beta,' persistence=',fit%margins(i)%alpha+fit%margins(i)%beta
   end do
   print '(a,f10.6)', 'DCC alpha: ', fit%dcc%spec%alpha(1)
   print '(a,f10.6)', 'DCC beta:  ', fit%dcc%spec%beta(1)
   print '(a,f14.4)', 'DCC log likelihood: ', fit%dcc%log_likelihood
   print '(a,i0)', 'status: ', fit%status

contains

   subroutine read_csv_with_header(path,x)
      character(len=*), intent(in) :: path
      real(dp), allocatable, intent(out) :: x(:,:)
      character(len=16384) :: line
      integer :: unit, ios, nrow, ncol, i

      open(newunit=unit,file=path,status='old',action='read',iostat=ios)
      if (ios /= 0) error stop 'could not open CSV file'
      read(unit,'(a)',iostat=ios) line
      if (ios /= 0) error stop 'CSV file is empty'
      ncol = 1+count_commas(trim(line))
      nrow = 0
      do
         read(unit,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) > 0) nrow = nrow+1
      end do
      rewind(unit)
      read(unit,'(a)') line
      allocate(x(nrow,ncol))
      do i = 1, nrow
         read(unit,'(a)',iostat=ios) line
         if (ios /= 0) error stop 'unexpected end of CSV file'
         call commas_to_spaces(line)
         read(line,*,iostat=ios) x(i,:)
         if (ios /= 0) error stop 'CSV must contain a header and numeric columns only'
      end do
      close(unit)
   end subroutine read_csv_with_header

   pure integer function count_commas(line) result(n)
      character(len=*), intent(in) :: line
      integer :: j
      n = 0
      do j = 1, len_trim(line)
         if (line(j:j) == ',') n = n+1
      end do
   end function count_commas

   subroutine commas_to_spaces(line)
      character(len=*), intent(inout) :: line
      integer :: j
      do j = 1, len_trim(line)
         if (line(j:j) == ',') line(j:j) = ' '
      end do
   end subroutine commas_to_spaces

end program fit_csv
