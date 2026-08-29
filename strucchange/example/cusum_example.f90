program cusum_example
   use r_kinds, only : dp
   use strucchange, only : ols_cusum, process_max
   implicit none
   integer, parameter :: n = 30
   real(dp) :: x(n, 2), y(n)
   real(dp), allocatable :: process(:, :)
   integer :: i, info

   do i = 1, n
      x(i, 1) = 1.0_dp
      x(i, 2) = real(i, dp)
      y(i) = 2.0_dp + 0.3_dp * real(i, dp) + 0.05_dp * cos(real(i, dp))
      if (i > 18) y(i) = y(i) + 0.8_dp
   end do

   call ols_cusum(x, y, process, info)
   if (info /= 0) error stop "OLS-CUSUM computation failed"
   print '(a,es12.4)', "Maximum absolute OLS-CUSUM: ", process_max(process)
end program cusum_example
