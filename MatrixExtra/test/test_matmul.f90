program test_matmul
   use matrixextra
   implicit none
   real(dp) :: x(4,3), y(3,2), d(2,4)
   real(dp), allocatable :: got(:,:), v(:)
   type(csr_matrix) :: a,b,c
   type(csc_matrix) :: cs
   integer :: info
   x=reshape([1.0_dp,0.0_dp,2.0_dp,0.0_dp, 0.0_dp,3.0_dp,0.0_dp,4.0_dp, 5.0_dp,0.0_dp,6.0_dp,0.0_dp],[4,3])
   y=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[3,2])
   call csr_from_dense(x,a)
   got=csr_dense_matmul(a,y)
   if (maxval(abs(got-matmul(x,y)))>1e-12_dp) error stop 1
   v=csr_matvec_extra(a,[2.0_dp,-1.0_dp,0.5_dp])
   if (maxval(abs(v-matmul(x,[2.0_dp,-1.0_dp,0.5_dp])))>1e-12_dp) error stop 2
   call csr_from_dense(y,b)
   call csr_csr_matmul(a,b,c,info)
   if (maxval(abs(csr_to_dense(c)-matmul(x,y)))>1e-12_dp) error stop 3
   d=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp,7.0_dp,8.0_dp],[2,4])
   call csc_from_csr(a,cs)
   got=dense_csc_matmul(d,cs)
   if (maxval(abs(got-matmul(d,x)))>1e-12_dp) error stop 4
   got=csr_crossprod(a,a)
   if (maxval(abs(got-matmul(transpose(x),x)))>1e-12_dp) error stop 5
   got=csr_tcrossprod(a,a)
   if (maxval(abs(got-matmul(x,transpose(x))))>1e-12_dp) error stop 6
   print *, 'PASS test_matmul'
end program
