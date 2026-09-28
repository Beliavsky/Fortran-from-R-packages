program basic_example
   use psych, only: dp, describe_vector, alpha_to_r, fisher_z, describe_result
   implicit none
   real(dp) :: x(5)
   type(describe_result) :: d

   x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   call describe_vector(x, d)
   write (*, '(a,f6.3)') 'mean = ', d%mean
   write (*, '(a,f6.3)') 'sd = ', d%sd
   write (*, '(a,f6.3)') 'alpha-to-r (alpha=.8, k=5) = ', alpha_to_r(0.8_dp, 5)
   write (*, '(a,f6.3)') 'Fisher z(.5) = ', fisher_z(0.5_dp)
end program basic_example
