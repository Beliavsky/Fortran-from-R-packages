program distribution_table
   use skellam, only : dp, i8, dskellam, pskellam
   implicit none
   integer(i8) :: k

   print '(a)', ' k          PMF             CDF'
   do k = -6_i8, 8_i8
      print '(i3,2(2x,f14.10))', k, dskellam(k, 3.0_dp, 2.0_dp), &
         pskellam(real(k, dp), 3.0_dp, 2.0_dp)
   end do
end program distribution_table
