! SPDX-License-Identifier: GPL-2.0-or-later
! Optional R validation driver, compiled in a temporary directory by run_reference.R.
program reference_driver
   use kendall_module, only : dp, kendall, mann_kendall, seasonal_mann_kendall, kendall_result
   implicit none
   character(len = 4096) :: input_file
   integer :: unit, operation, n, period, method, i
   real(dp), allocatable :: x(:), y(:)
   type(kendall_result) :: result

   call get_command_argument(1, input_file)
   open(newunit = unit, file = trim(input_file), status = 'old', action = 'read')
   read(unit, *) operation, n, period, method
   allocate(x(n), y(n))
   do i = 1, n
      read(unit, *) x(i), y(i)
   end do
   close(unit)
   select case (operation)
   case (1)
      call kendall(x, y, result, method)
   case (2)
      call mann_kendall(x, result, method)
   case (3)
      call seasonal_mann_kendall(x, period, result)
   case default
      error stop 'Unknown reference operation'
   end select
   write(*, '(3(i0,1x),5(es25.16e3,1x))') result%status, result%n_used, merge(1, 0, result%exact), &
      result%tau, result%sl, result%s, result%d, result%vars
end program reference_driver
