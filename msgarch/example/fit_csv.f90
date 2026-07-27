! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program fit_csv
   use msgarch
   implicit none
   character(len=256) :: filename, mode, model1, dist1, model2, dist2
   character(len=12) :: models(2)
   character(len=8) :: distributions(2)
   real(dp), allocatable :: y(:)
   type(msgarch_spec) :: spec
   type(fit_result) :: fit
   integer :: nargs, i

   nargs = command_argument_count()
   if (nargs < 1) then
      write(*,'(a)') 'usage: fit_csv FILE [single|markov|mixture] [MODEL1] [DIST1] [MODEL2] [DIST2]'
      error stop 1
   end if
   call get_command_argument(1, filename)
   mode = 'single'; model1 = 'sGARCH'; dist1 = 'norm'; model2 = 'gjrGARCH'; dist2 = 'std'
   if (nargs >= 2) call get_command_argument(2, mode)
   if (nargs >= 3) call get_command_argument(3, model1)
   if (nargs >= 4) call get_command_argument(4, dist1)
   if (nargs >= 5) call get_command_argument(5, model2)
   if (nargs >= 6) call get_command_argument(6, dist2)

   call read_series(trim(filename), y)
   models(1) = trim(model1); models(2) = trim(model2)
   distributions(1) = trim(dist1); distributions(2) = trim(dist2)
   select case (trim(mode))
   case ('single')
      spec = create_spec(models(1:1), distributions(1:1))
   case ('markov')
      spec = create_spec(models, distributions)
   case ('mixture')
      spec = create_spec(models, distributions, .true.)
   case default
      error stop 'fit_csv: mode must be single, markov, or mixture'
   end select

   fit = fit_ml(spec, y, max_iterations=500)
   write(*,'(a,i0)') 'observations: ', size(y)
   write(*,'(a,l1)') 'converged: ', fit%converged
   write(*,'(a,f16.6)') 'log likelihood: ', fit%loglik
   write(*,'(a,f16.6)') 'AIC: ', fit%aic
   write(*,'(a,f16.6)') 'BIC: ', fit%bic
   do i = 1, size(fit%parameters)
      write(*,'(i4,2x,es18.8,2x,es18.8)') i, fit%parameters(i), fit%standard_error(i)
   end do
contains
   subroutine read_series(path, values)
      character(len=*), intent(in) :: path
      real(dp), allocatable, intent(out) :: values(:)
      character(len=1024) :: line, token
      integer :: unit, ios, n, comma
      real(dp) :: value
      n = 0
      open(newunit=unit, file=path, status='old', action='read')
      do
         read(unit,'(a)',iostat=ios) line
         if (ios /= 0) exit
         comma = index(line, ',', back=.true.)
         if (comma > 0) then
            token = adjustl(line(comma+1:))
         else
            token = adjustl(line)
         end if
         read(token,*,iostat=ios) value
         if (ios == 0) n = n + 1
      end do
      rewind(unit)
      allocate(values(n)); n = 0
      do
         read(unit,'(a)',iostat=ios) line
         if (ios /= 0) exit
         comma = index(line, ',', back=.true.)
         if (comma > 0) then
            token = adjustl(line(comma+1:))
         else
            token = adjustl(line)
         end if
         read(token,*,iostat=ios) value
         if (ios == 0) then
            n = n + 1
            values(n) = value
         end if
      end do
      close(unit)
      if (size(values) < 20) error stop 'fit_csv: too few numeric observations'
   end subroutine read_series
end program fit_csv
