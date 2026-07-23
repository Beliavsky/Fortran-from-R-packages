! Part of the experimental modern Fortran translation of fGarch 4052.93.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original fGarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

program fit_csv
   use fgarch, only : dp, garch_fit_result, fit_garch11, dist_norm
   implicit none

   character(len=512) :: filename, line
   real(dp), allocatable :: x(:)
   real(dp) :: value
   integer :: unit, ios, n, i
   type(garch_fit_result) :: fit

   if (command_argument_count() < 1) then
      print '(a)', 'usage: fpm run --example fit_csv -- FILE.csv'
      stop 1
   end if
   call get_command_argument(1,filename)

   open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
   if (ios /= 0) error stop 'could not open input file'
   n = 0
   do
      read(unit,'(a)',iostat=ios) line
      if (ios /= 0) exit
      read(line,*,iostat=ios) value
      if (ios == 0) n = n+1
   end do
   if (n < 20) error stop 'need at least 20 numeric observations'
   rewind(unit)
   allocate(x(n))
   i = 0
   do
      read(unit,'(a)',iostat=ios) line
      if (ios /= 0) exit
      read(line,*,iostat=ios) value
      if (ios == 0) then
         i = i+1
         x(i) = value
      end if
   end do
   close(unit)

   fit = fit_garch11(x,cond_dist=dist_norm,max_iterations=1800)
   print '(a,i0)', 'observations: ', size(x)
   print '(a,f12.8)', 'mean:  ', fit%spec%mean
   print '(a,es14.6)', 'omega: ', fit%spec%omega
   print '(a,f10.6)', 'alpha: ', fit%spec%alpha(1)
   print '(a,f10.6)', 'beta:  ', fit%spec%beta(1)
   print '(a,f14.3)', 'log likelihood: ', fit%log_likelihood
   print '(a,a)', 'status: ', trim(fit%message)
end program fit_csv
