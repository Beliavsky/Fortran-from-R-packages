program test_sparse_matrix
   use trustoptim
   implicit none
   type(sparse_symmetric_matrix) :: h
   real(dp) :: a(3,3), b(3,3), x(3), y(3)

   a = reshape([4.0_dp,1.0_dp,0.0_dp, 1.0_dp,3.0_dp,2.0_dp, 0.0_dp,2.0_dp,5.0_dp],[3,3])
   call h%set_from_dense(a)
   if (h%nnz /= 5) error stop 'wrong lower-triangle nnz'
   x = [1.0_dp,2.0_dp,3.0_dp]
   call h%matvec(x,y)
   if (maxval(abs(y-matmul(a,x))) > 1.0e-12_dp) error stop 'sparse matvec mismatch'
   call h%to_dense(b)
   if (maxval(abs(a-b)) > 1.0e-12_dp) error stop 'sparse dense roundtrip mismatch'
   write(*,*) 'PASS sparse symmetric matrix'
end program test_sparse_matrix
