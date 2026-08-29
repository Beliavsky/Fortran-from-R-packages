program breakpoints_example
   use r_kinds, only : dp
   use strucchange, only : breakpoint_path_result
   use strucchange, only : best_break_count, compute_breakpoint_path
   implicit none
   integer, parameter :: n = 40
   real(dp) :: x(n, 2), y(n)
   type(breakpoint_path_result) :: path
   integer :: i, selected

   do i = 1, n
      x(i, 1) = 1.0_dp
      x(i, 2) = real(i, dp)
      if (i <= 20) then
         y(i) = 1.0_dp + 0.10_dp * real(i, dp) + 0.04_dp * sin(real(i, dp))
      else
         y(i) = -1.0_dp + 0.24_dp * real(i, dp) + 0.04_dp * sin(real(i, dp))
      end if
   end do

   call compute_breakpoint_path(x, y, 8, 3, path)
   if (path%info /= 0) error stop "breakpoint computation failed"
   selected = best_break_count(path)
   print '(a,i0)', "BIC-selected number of breaks: ", selected
   if (selected > 0) print '(a,*(i0,1x))', "Breakpoints: ", path%breakpoints(1:selected, selected)
end program breakpoints_example
