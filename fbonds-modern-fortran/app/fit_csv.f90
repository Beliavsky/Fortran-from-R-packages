! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! Based on fBonds, copyright its original authors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License, version 2 or later.
program fit_csv
   use fbonds_kinds, only : dp
   use fbonds_csv, only : read_maturity_rate_csv
   use fbonds_term_structure, only : term_structure_fit, fit_nelson_siegel, fit_svensson
   implicit none
   character(len=1024) :: filename
   character(len=32) :: model, objective
   real(dp), allocatable :: maturity(:), rate(:)
   type(term_structure_fit) :: fit = term_structure_fit()
   integer :: status, i

   if (command_argument_count() < 1) then
      print '(a)', 'usage: fit_csv FILE [ns|svensson] [l1|sse]'
      error stop 2
   end if
   call get_command_argument(1, filename)
   model = 'ns'
   objective = 'l1'
   if (command_argument_count() >= 2) call get_command_argument(2, model)
   if (command_argument_count() >= 3) call get_command_argument(3, objective)

   call read_maturity_rate_csv(trim(filename), maturity, rate, status)
   if (status /= 0 .or. size(rate) < 4) then
      print '(a)', 'could not read at least four maturity,rate rows'
      error stop 3
   end if

   select case (trim(lower_string(model)))
   case ('ns', 'nelson-siegel', 'nelson_siegel')
      call fit_nelson_siegel(rate, maturity, fit)
   case ('svensson', 'nss')
      call fit_svensson(rate, maturity, fit, objective=trim(objective))
   case default
      print '(a)', 'model must be ns or svensson'
      error stop 4
   end select

   print '(a,a)', 'model: ', trim(fit%model)
   print '(a,a)', 'objective: ', trim(fit%objective)
   print '(a,l1)', 'converged: ', fit%converged
   print '(a,i0)', 'iterations: ', fit%iterations
   print '(a,i0)', 'evaluations: ', fit%evaluations
   print '(a,es18.10)', 'objective_value: ', fit%objective_value
   print '(a,es18.10)', 'sse: ', fit%sse
   print '(a,es18.10)', 'mae: ', fit%mae
   print '(a,es18.10)', 'rmse: ', fit%rmse
   do i = 1, size(fit%parameters)
      print '(a,i0,a,es18.10)', 'parameter(', i, '): ', fit%parameters(i)
   end do
contains
   pure function lower_string(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code
      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lower_string
end program fit_csv
