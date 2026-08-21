program test_ops_utils
   use matrixextra
   implicit none
   type(csr_matrix) :: a,b,c
   real(dp) :: x(3,3), y(3,3), want(3,3)
   integer :: info
   x=reshape([1.0_dp,0.0_dp,2.0_dp,0.0_dp,-3.0_dp,0.0_dp,4.0_dp,0.0_dp,5.0_dp],[3,3])
   y=reshape([0.0_dp,6.0_dp,2.0_dp,0.0_dp,3.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp],[3,3])
   call csr_from_dense(x,a); call csr_from_dense(y,b)
   call csr_elem_add(a,b,c,info)
   if (maxval(abs(csr_to_dense(c)-(x+y)))>1e-13_dp) error stop 1
   call csr_elem_subtract(a,b,c,info)
   if (maxval(abs(csr_to_dense(c)-(x-y)))>1e-13_dp) error stop 2
   call csr_elem_multiply(a,b,c,info)
   if (maxval(abs(csr_to_dense(c)-(x*y)))>1e-13_dp) error stop 3
   call csr_logical_and(a,b,c,info)
   want=0.0_dp; where(abs(x)>0.0_dp .and. abs(y)>0.0_dp) want=1.0_dp
   if (maxval(abs(csr_to_dense(c)-want))>1e-13_dp) error stop 4
   call csr_logical_or(a,b,c,info)
   want=0.0_dp; where(abs(x)>0.0_dp .or. abs(y)>0.0_dp) want=1.0_dp
   if (maxval(abs(csr_to_dense(c)-want))>1e-13_dp) error stop 5
   call csr_from_dense(x,a)
   call csr_multiply_vector(a,[2.0_dp,3.0_dp])
   want=x
   want(1,1)=x(1,1)*2.0_dp; want(3,1)=x(3,1)*2.0_dp
   want(2,2)=x(2,2)*2.0_dp; want(1,3)=x(1,3)*2.0_dp; want(3,3)=x(3,3)*2.0_dp
   ! R column-major recycling: flattened positions 1,3,5,7,9 are all odd -> multiplier 2.
   if (maxval(abs(csr_to_dense(a)-want))>1e-13_dp) error stop 6
   call csr_from_dense(x,a); call csr_apply_abs(a)
   if (maxval(abs(csr_to_dense(a)-abs(x)))>1e-13_dp) error stop 7
   call csr_remove_zeros(a)
   if (.not. csr_check(a)) error stop 8
   print *, 'PASS test_ops_utils'
end program
