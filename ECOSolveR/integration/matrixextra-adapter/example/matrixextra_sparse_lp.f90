program matrixextra_sparse_lp
    use matrix_kinds, only : mdp => dp
    use matrix_sparse, only : csr_matrix, csc_matrix, csr_from_triplet
    use matrixextra_conversions, only : csr_to_csc
    use ecos, only : dp, ecos_dims, ecos_problem, ecos_result, ecos_csolve, ECOS_OPTIMAL
    use ecos_matrixextra_adapter, only : setup_problem_matrixextra
    implicit none
    type(csr_matrix) :: gr
    type(csc_matrix) :: gc
    type(ecos_dims) :: dims
    type(ecos_problem) :: prob
    type(ecos_result) :: result
    integer :: rows(3),cols(3),info
    real(mdp) :: vals(3)
    real(dp) :: c(3),h(3)

    rows=[1,2,3]; cols=[1,2,3]; vals=-1.0_mdp
    call csr_from_triplet(3,3,rows,cols,vals,gr,info)
    if(info/=0) error stop 'Matrix sparse construction failed'
    call csr_to_csc(gr,gc)
    c=1.0_dp; h=-1.0_dp; dims%l=3
    call setup_problem_matrixextra(prob,c,gc,h,dims,ierr=info)
    if(info/=0) error stop 'ECOS setup failed'
    call ecos_csolve(prob,result)
    if(result%exitflag/=ECOS_OPTIMAL) error stop 'ECOS solve failed'
    if(maxval(abs(result%x-1.0_dp))>1.0e-7_dp) error stop 'wrong sparse LP solution'
    if(allocated(prob%gmat)) error stop 'adapter unexpectedly densified G'
    print '(a,3f10.5)', 'MatrixExtra sparse LP x =',result%x
end program matrixextra_sparse_lp
