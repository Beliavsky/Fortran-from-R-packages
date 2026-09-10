! SPDX-License-Identifier: GPL-2.0-only
program named_constraints
   use marss_api
   implicit none

   real(dp) :: fixed_values(3)
   real(dp) :: starts(2)
   character(len=64) :: expressions(3)
   logical :: is_expression(3)
   type(marss_constraint_block) :: block
   type(marss_constraints) :: constraints
   character(len=64), allocatable :: names(:)
   integer :: info
   integer :: i

   fixed_values = [5.0_dp, 0.0_dp, 0.0_dp]
   expressions = [character(len=64) :: "", "alpha", "2+0.5*alpha+3*beta"]
   is_expression = [.false., .true., .true.]
   starts = [0.25_dp, -0.5_dp]

   call marss_constraint_from_entries(fixed_values, expressions, is_expression, block, info, starts)
   if (info /= 0) error stop "constraint expression conversion failed"

   constraints%a = block
   call marss_vectorized_parameter_names(constraints, names)
   do i = 1, size(names)
      write (*, '(a)') trim(names(i))
   end do
end program named_constraints
