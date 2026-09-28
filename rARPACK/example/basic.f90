program basic
   use rarpack, only: dp, eigs_sym, eigs_sym_result, svds, svds_result
   implicit none

   real(dp) :: a(8, 8)
   real(dp) :: b(6, 4)
   type(eigs_sym_result) :: er
   type(svds_result) :: sr
   integer :: i

   a = 0.0_dp
   do i = 1, 8
      a(i, i) = real(i, dp)
   end do
   er = eigs_sym(a, 3, which = 'LA')
   if (er%info /= 0) error stop 'eigs_sym failed'
   print '(a,3f10.4)', 'largest eigenvalues: ', er%values

   b = 0.0_dp
   b(1, 1) = 4.0_dp
   b(2, 2) = 3.0_dp
   b(3, 3) = 2.0_dp
   b(4, 4) = 1.0_dp
   sr = svds(b, 2)
   if (sr%info /= 0) error stop 'svds failed'
   print '(a,2f10.4)', 'largest singular values: ', sr%d
end program basic
