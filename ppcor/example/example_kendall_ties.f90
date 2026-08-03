program example_kendall_ties
   use ppcor, only : dp, ppcor_result, pcor, ppcor_kendall
   implicit none
   real(dp) :: x(12,4)
   type(ppcor_result) :: result

   x(:,1) = real([1,1,2,2,3,3,4,5,5,6,7,8],dp)
   x(:,2) = real([4,3,3,2,2,1,1,2,4,5,6,6],dp)
   x(:,3) = real([2,1,2,3,3,4,4,5,6,5,7,8],dp)
   x(:,4) = real([1,2,2,3,4,4,5,5,6,7,7,8],dp)

   call pcor(x, result, ppcor_kendall)
   print '(a,f10.6)', 'Kendall tau-b partial correlation = ', result%estimate(1,2)
   print '(a,f10.6)', 'normal statistic = ', result%statistic(1,2)
   print '(a,es11.4)', 'two-sided p-value = ', result%p_value(1,2)
end program example_kendall_ties
