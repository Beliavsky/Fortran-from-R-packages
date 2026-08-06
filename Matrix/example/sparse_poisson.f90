! SPDX-License-Identifier: GPL-3.0-only
program sparse_poisson
   use matrix, only : dp, csr_matrix, csr_from_triplet, csr_matvec, matrix_success
   implicit none
   type(csr_matrix) :: a
   integer, allocatable :: rows(:), cols(:)
   real(dp), allocatable :: values(:), x(:), y(:)
   integer :: n, i, k, info

   n = 5
   allocate(rows(3 * n - 2), cols(3 * n - 2), values(3 * n - 2))
   k = 0
   do i = 1, n
      if (i > 1) then
         k = k + 1
         rows(k) = i
         cols(k) = i - 1
         values(k) = -1.0_dp
      end if
      k = k + 1
      rows(k) = i
      cols(k) = i
      values(k) = 2.0_dp
      if (i < n) then
         k = k + 1
         rows(k) = i
         cols(k) = i + 1
         values(k) = -1.0_dp
      end if
   end do
   call csr_from_triplet(n, n, rows, cols, values, a, info)
   if (info /= matrix_success) error stop 'construction failed'
   x = [(real(i, dp), i = 1, n)]
   call csr_matvec(a, x, y, info)
   print '(a,i0)', 'nonzeros: ', a%nnz()
   print '(a,*(1x,f7.2))', 'A*x:', y
end program sparse_poisson
