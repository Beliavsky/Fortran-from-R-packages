program test_conversions_slice
   use matrixextra
   implicit none
   real(dp) :: x(4,5), y(3,3)
   type(csr_matrix) :: a,b
   type(csc_matrix) :: c
   type(coo_matrix) :: t
   type(sparse_vector) :: v
   integer :: rows(3),cols(3), info
   x=reshape([1.0_dp,0.0_dp,4.0_dp,0.0_dp, &
              0.0_dp,2.0_dp,0.0_dp,7.0_dp, &
              3.0_dp,0.0_dp,5.0_dp,0.0_dp, &
              0.0_dp,0.0_dp,6.0_dp,8.0_dp, &
              9.0_dp,0.0_dp,0.0_dp,10.0_dp],[4,5])
   call csr_from_dense(x,a)
   if (.not. a%valid()) error stop 1
   call coo_from_csr(a,t)
   call csr_from_coo(t,b,info)
   if (maxval(abs(csr_to_dense(b)-x))>1e-13_dp) error stop 2
   call csc_from_coo(t,c,info)
   call csc_to_csr(c,b)
   if (maxval(abs(csr_to_dense(b)-x))>1e-13_dp) error stop 3
   rows=[4,2,2]; cols=[5,1,3]
   call csr_slice(a,rows,cols,b)
   y=reshape([x(4,5),x(2,5),x(2,5),x(4,1),x(2,1),x(2,1), &
      x(4,3),x(2,3),x(2,3)],[3,3])
   ! reshape above is column major and row/col selectors define B(i,j)=X(rows(i),cols(j))
   y=0.0_dp
   y(1,:)=[x(4,5),x(4,1),x(4,3)]
   y(2,:)=[x(2,5),x(2,1),x(2,3)]
   y(3,:)=[x(2,5),x(2,1),x(2,3)]
   if (maxval(abs(csr_to_dense(b)-y))>1e-13_dp) error stop 4
   call sparse_vector_from_dense([0.0_dp,2.0_dp,0.0_dp,-3.0_dp],v)
   if (v%nnz()/=2 .or. maxval(abs(sparse_vector_to_dense(v)- &
      [0.0_dp,2.0_dp,0.0_dp,-3.0_dp]))>1e-13_dp) error stop 5
   print *, 'PASS test_conversions_slice'
end program
