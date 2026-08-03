program sparse_triplet_example
   use highs
   implicit none
   type(highs_sparse_matrix) :: matrix
   real(dp), allocatable :: dense(:,:)

   matrix = highs_csc_from_triplets(3, 3, [1,3,2,1,1], [1,1,2,3,3], &
      [2.0_dp,4.0_dp,3.0_dp,0.5_dp,1.5_dp])
   dense = matrix%to_dense()
   print '(a,i0)', "nonzeros after duplicate aggregation: ", matrix%nnz()
   print '(3(f8.2,1x))', transpose(dense)
end program sparse_triplet_example
