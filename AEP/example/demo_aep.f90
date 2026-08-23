program demo_aep
   use aep
   implicit none
   real(dp) :: x(10)
   call raep(x, alpha=1.5_dp, sigma=1.0_dp, mu=0.0_dp, epsilon=0.2_dp)
   print '(10f9.4)', x
end program
