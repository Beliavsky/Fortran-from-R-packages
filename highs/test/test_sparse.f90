program test_sparse
   use highs
   implicit none
   real(dp) :: a(3,3), b(3,3)
   type(highs_sparse_matrix) :: csc, csr, tri
   integer :: row(4), col(4)
   real(dp) :: val(4)

   a = reshape([1.0_dp,0.0_dp,4.0_dp, 0.0_dp,3.0_dp,0.0_dp, 2.0_dp,0.0_dp,5.0_dp], [3,3])
   csc = highs_csc_from_dense(a)
   csr = highs_csr_from_dense(a)
   if (.not. csc%valid() .or. .not. csr%valid()) error stop "invalid sparse matrix"
   if (csc%nnz() /= 5 .or. csr%nnz() /= 5) error stop "wrong nnz"
   b = csc%to_dense()
   if (maxval(abs(a-b)) > 0.0_dp) error stop "CSC round trip failed"
   b = csr%to_dense()
   if (maxval(abs(a-b)) > 0.0_dp) error stop "CSR round trip failed"

   row = [1,1,3,3]
   col = [1,1,2,3]
   val = [0.25_dp,0.75_dp,4.0_dp,5.0_dp]
   tri = highs_csc_from_triplets(3,3,row,col,val)
   b = tri%to_dense()
   if (maxval(abs([b(1,1)-1.0_dp, b(3,2)-4.0_dp, b(3,3)-5.0_dp])) > 0.0_dp) &
      error stop "triplet aggregation failed"
   print *, "test_sparse: PASS"
end program test_sparse
