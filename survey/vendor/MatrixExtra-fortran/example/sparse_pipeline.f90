program sparse_pipeline
   use matrixextra
   implicit none
   type(csr_matrix) :: a,b
   real(dp) :: x(4,4)
   integer :: rows(3),cols(2)
   x=reshape([1.0_dp,0.0_dp,4.0_dp,0.0_dp, &
              0.0_dp,2.0_dp,0.0_dp,0.0_dp, &
              3.0_dp,0.0_dp,5.0_dp,0.0_dp, &
              0.0_dp,6.0_dp,0.0_dp,7.0_dp],[4,4])
   call csr_from_dense(x,a)
   rows=[4,1,3]; cols=[4,1]
   call csr_slice(a,rows,cols,b)
   print '(a,i0)', 'nnz = ', b%nnz()
   print '(a,es14.6)', 'Frobenius norm = ', csr_norm(b,'F')
   print *, csr_to_dense(b)
end program
