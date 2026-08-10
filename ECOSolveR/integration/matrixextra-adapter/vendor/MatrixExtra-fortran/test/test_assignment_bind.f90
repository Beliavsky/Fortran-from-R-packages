program test_assignment_bind
   use matrixextra
   implicit none
   real(dp) :: x(3,3)
   type(csr_matrix) :: a,b,c
   integer :: info
   x=reshape([1.0_dp,0.0_dp,4.0_dp,0.0_dp,2.0_dp,0.0_dp,3.0_dp,0.0_dp,5.0_dp],[3,3])
   call csr_from_dense(x,a)
   call csr_set_value(a,2,1,7.0_dp); x(2,1)=7.0_dp
   call csr_set_row_constant(a,1,2.0_dp); x(1,:)=2.0_dp
   call csr_set_col_constant(a,3,0.0_dp); x(:,3)=0.0_dp
   if (maxval(abs(csr_to_dense(a)-x))>1e-13_dp) error stop 1
   call csr_from_dense(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2]),a)
   call csr_from_dense(reshape([5.0_dp,6.0_dp],[1,2]),b)
   call csr_rbind(a,b,c,info)
   if (maxval(abs(csr_to_dense(c)-reshape([1.0_dp,2.0_dp,5.0_dp,3.0_dp,4.0_dp,6.0_dp],[3,2])))>1e-13_dp) error stop 2
   call csr_from_dense(reshape([1.0_dp,2.0_dp],[2,1]),a)
   call csr_from_dense(reshape([3.0_dp,4.0_dp,5.0_dp,6.0_dp],[2,2]),b)
   call csr_cbind(a,b,c,info)
   if (maxval(abs(csr_to_dense(c)-reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[2,3])))>1e-13_dp) error stop 3
   print *, 'PASS test_assignment_bind'
end program
