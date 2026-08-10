program matmul_example
   use matrixextra
   implicit none
   type(csr_matrix) :: a
   real(dp) :: x(3,3), b(3,2)
   real(dp), allocatable :: c(:,:)
   x=reshape([1.0_dp,0.0_dp,2.0_dp,0.0_dp,3.0_dp,0.0_dp,4.0_dp,0.0_dp,5.0_dp],[3,3])
   b=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[3,2])
   call csr_from_dense(x,a)
   c=csr_dense_matmul(a,b)
   print *, c
end program
