! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program fit_csv
   use gogarch
   implicit none
   character(len=1024) :: filename, method, line
   character(len=128) :: date_token, argument
   real(dp), allocatable :: data(:,:), mean_forecast(:,:), covariance_forecast(:,:,:)
   real(dp), allocatable :: means(:), omegas(:), arch(:,:), leverage(:,:), garch(:,:), delta(:), shape(:), skew(:)
   type(gogarch_fit) :: fit
   type(univariate_spec) :: spec
   integer :: unit, ios, n, m, i, j, k

   if (command_argument_count() < 1) then
      write(*,'(a)') 'usage: fit_csv FILE [ica|mm|nls|ml] [garch|aparch] [distribution] [p] [o] [q]'
      write(*,'(a)') 'distributions: norm snorm std sstd ged sged'
      write(*,'(a)') 'CSV format: Date,Asset1,Asset2,... with one header row.'
      error stop 2
   end if
   call get_command_argument(1,filename)
   method = 'ica'
   if (command_argument_count() >= 2) call get_command_argument(2,method)
   spec = univariate_spec()
   if (command_argument_count() >= 3) call get_command_argument(3,spec%model)
   if (command_argument_count() >= 4) call get_command_argument(4,spec%distribution)
   if (command_argument_count() >= 5) then
      call get_command_argument(5,argument)
      read(argument,*,iostat=ios) spec%p
      if (ios /= 0) error stop 'invalid ARCH order p'
   end if
   if (command_argument_count() >= 6) then
      call get_command_argument(6,argument)
      read(argument,*,iostat=ios) spec%o
      if (ios /= 0) error stop 'invalid leverage order o'
   end if
   if (command_argument_count() >= 7) then
      call get_command_argument(7,argument)
      read(argument,*,iostat=ios) spec%q
      if (ios /= 0) error stop 'invalid GARCH order q'
   end if
   if (trim(adjustl(spec%model)) == 'aparch') then
      spec%delta = 2.0_dp
      spec%fit_delta = .true.
   end if
   select case (trim(adjustl(spec%distribution)))
   case ('std','sstd')
      spec%shape = 8.0_dp
   case ('ged','sged')
      spec%shape = 1.5_dp
   end select
   if (trim(adjustl(spec%distribution)) == 'snorm' .or. trim(adjustl(spec%distribution)) == 'sstd' .or. &
       trim(adjustl(spec%distribution)) == 'sged') spec%skew = 1.1_dp
   if (.not. validate_specification(spec)) error stop 'invalid univariate model specification'

   open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
   if (ios /= 0) error stop 'could not open CSV file'
   read(unit,'(a)',iostat=ios) line
   if (ios /= 0) error stop 'could not read CSV header'
   m = count_commas(trim(line))
   if (m < 1) error stop 'CSV must contain a date column and at least one asset column'
   n = 0
   do
      read(unit,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) n = n+1
   end do
   if (n < 10) error stop 'CSV must contain at least ten observations'
   rewind(unit)
   read(unit,'(a)') line
   allocate(data(n,m))
   i = 0
   do
      read(unit,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      i = i+1
      call commas_to_spaces(line)
      read(line,*,iostat=ios) date_token,data(i,:)
      if (ios /= 0) error stop 'could not parse a CSV data row'
   end do
   close(unit)

   fit = fit_gogarch(data,trim(adjustl(method)),lag_max=2,max_outer_iterations=60,max_garch_iterations=350, &
      factor_spec=spec)
   if (fit%status > 1) then
      write(*,'(a,i0)') 'fit failed with status ',fit%status
      error stop 3
   end if
   allocate(means(m),omegas(m),arch(m,spec%p),leverage(m,spec%o),garch(m,spec%q))
   allocate(delta(m),shape(m),skew(m),mean_forecast(5,m),covariance_forecast(m,m,5))
   call factor_coefficients_full(fit,means,omegas,arch,leverage,garch,delta,shape,skew)
   call forecast_gogarch(fit,5,mean_forecast,covariance_forecast)

   write(*,'(a,a)') 'method: ',trim(fit%method)
   write(*,'(a,a)') 'factor model: ',trim(fit%factor_spec%model)
   write(*,'(a,a)') 'conditional distribution: ',trim(fit%factor_spec%distribution)
   write(*,'(a,3(i0,1x))') 'orders p o q: ',spec%p,spec%o,spec%q
   write(*,'(a,i0)') 'observations: ',n
   write(*,'(a,i0)') 'assets: ',m
   write(*,'(a,f16.6)') 'log likelihood: ',fit%log_likelihood
   write(*,'(a,es12.4)') 'maximum reconstruction error: ',reconstruction_error(fit)
   do j = 1, m
      write(*,'(a,i0)') 'factor ',j
      write(*,'(a,1x,es16.8)') '  mean:',means(j)
      write(*,'(a,1x,es16.8)') '  omega:',omegas(j)
      write(*,'(a)',advance='no') '  arch:'
      do k = 1, spec%p
         write(*,'(1x,es16.8)',advance='no') arch(j,k)
      end do
      write(*,*)
      if (spec%o > 0) then
         write(*,'(a)',advance='no') '  leverage:'
         do k = 1, spec%o
            write(*,'(1x,es16.8)',advance='no') leverage(j,k)
         end do
         write(*,*)
      end if
      if (spec%q > 0) then
         write(*,'(a)',advance='no') '  garch:'
         do k = 1, spec%q
            write(*,'(1x,es16.8)',advance='no') garch(j,k)
         end do
         write(*,*)
      end if
      write(*,'(a,3(1x,es16.8))') '  delta shape skew:',delta(j),shape(j),skew(j)
   end do
   write(*,'(a)') 'one-step covariance forecast:'
   do j = 1, m
      write(*,'(*(1x,f13.7))') covariance_forecast(j,:,1)
   end do

contains

   pure integer function count_commas(text) result(number)
      character(len=*), intent(in) :: text
      integer :: k
      number = 0
      do k = 1, len_trim(text)
         if (text(k:k) == ',') number = number+1
      end do
   end function count_commas

   pure subroutine commas_to_spaces(text)
      character(len=*), intent(inout) :: text
      integer :: k
      do k = 1, len(text)
         if (text(k:k) == ',') text(k:k) = ' '
      end do
   end subroutine commas_to_spaces

end program fit_csv
