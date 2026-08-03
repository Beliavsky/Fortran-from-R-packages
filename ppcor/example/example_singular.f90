program example_singular
   use ppcor, only : dp, ppcor_result, pcor, ppcor_pearson
   implicit none
   real(dp) :: x(18,4)
   type(ppcor_result) :: result
   integer :: i
   real(dp) :: t

   do i = 1, size(x,1)
      t = real(i,dp)
      x(i,1) = sin(0.3_dp*t)
      x(i,2) = cos(0.2_dp*t)
      x(i,3) = x(i,1)
      x(i,4) = sin(0.7_dp*t) + 0.1_dp*t
   end do

   call pcor(x, result, ppcor_pearson)
   print '(a,l1)', 'used Moore-Penrose pseudoinverse: ', result%used_pseudoinverse
   print '(a,i0)', 'association-matrix rank: ', result%rank
   print '(a)', trim(result%message)
end program example_singular
