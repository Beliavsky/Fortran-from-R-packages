! Part of the experimental modern Fortran translation of tseries 0.10-62.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original tseries authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

program fit_csv
   use tseries, only : dp, test_result, arma_result, garch_result, jarque_bera_test, arma_fit, garch_fit
   implicit none
   character(len=1024) :: filename,line
   real(dp), allocatable :: x(:)
   real(dp) :: value
   integer :: unit,ios,n,i
   type(arma_result) :: ar
   type(garch_result) :: ga
   type(test_result) :: jb

   if(command_argument_count()<1) then
      print '(a)','usage: fit_csv FILE.csv'
      stop 1
   end if
   call get_command_argument(1,filename)
   open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
   if(ios/=0) error stop 'could not open input file'
   n=0
   do
      read(unit,'(a)',iostat=ios) line
      if(ios/=0) exit
      read(line,*,iostat=ios) value
      if(ios==0) n=n+1
   end do
   rewind(unit)
   allocate(x(n)); i=0
   do
      read(unit,'(a)',iostat=ios) line
      if(ios/=0) exit
      read(line,*,iostat=ios) value
      if(ios==0) then
         i=i+1; x(i)=value
      end if
   end do
   close(unit)
   if(n<10) error stop 'input needs at least ten numeric rows'

   ar=arma_fit(x,1,1,max_iterations=1000)
   ga=garch_fit(x-sum(x)/real(n,dp),1,1,max_iterations=1000)
   jb=jarque_bera_test(x)
   print '(a,i0)','observations: ',n
   print '(a,f12.6)','Jarque-Bera: ',jb%statistic
   print '(a,*(f12.6,1x))','ARMA(1,1): ',ar%coefficients
   print '(a,*(f12.6,1x))','GARCH(1,1): ',ga%coefficients
end program fit_csv
