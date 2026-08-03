! SPDX-License-Identifier: MIT
program compounding_example
   use fixedincome
   implicit none
   real(dp) :: rate, years

   rate = 0.05_dp
   years = 2.0_dp
   print '(a,f12.8)', 'simple:     ', compound(COMPOUND_SIMPLE, years, rate)
   print '(a,f12.8)', 'discrete:   ', compound(COMPOUND_DISCRETE, years, rate)
   print '(a,f12.8)', 'continuous: ', compound(COMPOUND_CONTINUOUS, years, rate)
end program compounding_example
