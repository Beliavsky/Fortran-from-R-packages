program test_linalg_transpose
   use matrixextra
   implicit none
   real(dp) :: x(3,4), d(3), f
   type(csr_matrix) :: a,b
   type(csc_matrix) :: c
   x=reshape([1.0_dp,0.0_dp,2.0_dp, 0.0_dp,3.0_dp,0.0_dp, 4.0_dp,0.0_dp,5.0_dp, 0.0_dp,6.0_dp,0.0_dp],[3,4])
   call csr_from_dense(x,a)
   if (abs(csr_norm(a,'1')-maxval(sum(abs(x),dim=1)))>1e-12_dp) error stop 1
   if (abs(csr_norm(a,'I')-maxval(sum(abs(x),dim=2)))>1e-12_dp) error stop 2
   if (abs(csr_norm(a,'F')-sqrt(sum(x*x)))>1e-12_dp) error stop 3
   if (abs(csr_norm(a,'M')-maxval(abs(x)))>1e-12_dp) error stop 4
   d=csr_diag(a)
   if (maxval(abs(d-[1.0_dp,3.0_dp,5.0_dp]))>1e-12_dp) error stop 5
   call csr_set_diag(a,[7.0_dp,8.0_dp,9.0_dp]); x(1,1)=7.0_dp; x(2,2)=8.0_dp; x(3,3)=9.0_dp
   if (maxval(abs(csr_to_dense(a)-x))>1e-12_dp) error stop 6
   call csr_transpose_shallow(a,c)
   call csc_transpose_shallow(c,b)
   if (maxval(abs(csr_to_dense(b)-x))>1e-12_dp) error stop 7
   f=csr_norm(a,'2')
   if (f<=0.0_dp) error stop 8
   print *, 'PASS test_linalg_transpose'
end program
