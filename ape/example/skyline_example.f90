program skyline_example
   use ape, only : dp, skyline_from_intervals, skyline_result
   implicit none

   integer :: lineages(4)
   real(dp) :: widths(4)
   type(skyline_result) :: result
   integer :: info
   integer :: i

   lineages = [5, 4, 3, 2]
   widths = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp]
   call skyline_from_intervals(lineages, widths, 0.35_dp, result, info)
   if (info /= 0) error stop 'skyline calculation failed'

   print '(a)', 'time population_size'
   do i = 1, result%parameter_count
      print '(f9.4,1x,f12.6)', result%time(i), result%population_size(i)
   end do
   print '(a,f12.6)', 'log likelihood: ', result%log_likelihood
   print '(a,f12.6)', 'AICc-corrected log likelihood: ', result%log_likelihood_aicc
end program skyline_example
