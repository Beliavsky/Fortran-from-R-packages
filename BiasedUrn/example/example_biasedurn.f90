program example_biasedurn
   use biasedurn
   implicit none
   integer, allocatable :: draws(:)
   integer :: x

   print '(a)', 'Fisher and Wallenius noncentral hypergeometric PMFs'
   print '(a)', ' x        Fisher              Wallenius'
   do x = 0, 8
      print '(i2,2(2x,f18.12))', x, &
         dfnchypergeo(x, 10, 15, 8, 2.5_dp), &
         dwnchypergeo(x, 10, 15, 8, 2.5_dp)
   end do

   call biasedurn_seed(2026)
   draws = rwnchypergeo(10, 10, 15, 8, 2.5_dp)
   print '(a,10(1x,i0))', 'Wallenius draws:', draws
end program example_biasedurn
