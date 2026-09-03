program birthdeath_example
   use ape, only : birthdeath_from_times, birthdeath_result, dp
   implicit none

   real(dp) :: times(4)
   type(birthdeath_result) :: fit
   integer :: info

   times = [3.613551448368_dp, 3.246318433096_dp, 1.406539523638_dp, 0.085104551586_dp]
   call birthdeath_from_times(times, fit, info)
   if (info /= 0) error stop 'birth-death fit failed'

   print '(a,f12.7)', 'd/b: ', fit%death_birth_ratio
   print '(a,f12.7)', 'b-d: ', fit%net_diversification
   print '(a,f12.7)', 'deviance: ', fit%deviance
   print '(a,2(1x,f12.7))', 'd/b 95% interval:', fit%confidence_interval(1, :)
   print '(a,2(1x,f12.7))', 'b-d 95% interval:', fit%confidence_interval(2, :)
end program birthdeath_example
