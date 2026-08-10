program basic
    use slam
    implicit none
    real(dp) :: a(3,3)
    real(dp), allocatable :: rs(:), gram(:,:)
    type(simple_triplet_matrix) :: x

    a=reshape([1.0_dp,0.0_dp,2.0_dp, &
               0.0_dp,3.0_dp,0.0_dp, &
               4.0_dp,0.0_dp,5.0_dp],[3,3])
    x=dense_to_stm(a)
    rs=x%row_sums()
    gram=stm_crossprod(x)

    print '(a,i0)', 'nonzeros: ',x%nnz()
    print '(a,3f8.3)', 'row sums: ',rs
    print '(a)', 'crossprod:'
    print '(3f10.3)',gram
end program basic
