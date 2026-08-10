! SPDX-License-Identifier: GPL-3.0-only
program sparse_ordering
   use matrix, only : dp, csr_matrix, csr_from_dense, reverse_cuthill_mckee, &
      minimum_degree_ordering, matrix_success
   implicit none
   type(csr_matrix) :: a
   real(dp), allocatable :: dense(:,:)
   integer, allocatable :: rcm(:), md(:)
   integer :: info

   allocate(dense(6, 6), source=0.0_dp)
   dense(1, [1, 4]) = 1.0_dp
   dense(2, [2, 5]) = 1.0_dp
   dense(3, [3, 6]) = 1.0_dp
   dense(4, [1, 4, 5]) = 1.0_dp
   dense(5, [2, 4, 5, 6]) = 1.0_dp
   dense(6, [3, 5, 6]) = 1.0_dp
   dense = max(dense, transpose(dense))
   call csr_from_dense(dense, a)
   call reverse_cuthill_mckee(a, rcm, info)
   if (info /= matrix_success) error stop 'RCM failed'
   call minimum_degree_ordering(a, md, info)
   print '(a,*(1x,i0))', 'RCM:', rcm
   print '(a,*(1x,i0))', 'minimum degree:', md
end program sparse_ordering
