program sparse_large_lp
    use ecos
    implicit none
    integer, parameter :: n=10000
    type(ecos_dims) :: dims
    type(ecos_csc_matrix) :: g
    type(ecos_result) :: result
    real(dp), allocatable :: c(:),h(:)
    integer :: i

    allocate(c(n),h(n),g%colptr(n+1),g%rowind(n),g%values(n))
    c=1.0_dp; h=-1.0_dp; dims%l=n
    g%nrow=n; g%ncol=n
    do i=1,n
        g%colptr(i)=i
        g%rowind(i)=i
        g%values(i)=-1.0_dp
    end do
    g%colptr(n+1)=n+1

    call ecos_csolve(c,g,h,dims,result)
    if(result%exitflag/=ECOS_OPTIMAL) error stop 'sparse solve failed'
    print '(a,i0)', 'variables: ',n
    print '(a,es14.6)', 'max |x-1|: ',maxval(abs(result%x-1.0_dp))
    print '(a,i0)', 'KKT upper nnz: ',result%kkt_nnz
    print '(a,i0)', 'LDL off-diagonal nnz: ',result%ldl_nnz
    print '(a,i0)', 'symbolic analyses: ',result%symbolic_analyses
    print '(a,i0)', 'numeric factorizations: ',result%numeric_factorizations
    print '(a,f8.3)', 'LDL/KKT fill ratio: ',result%factor_fill_ratio
    print '(a,es12.4)', 'ordering CPU seconds: ',result%time_ordering
    print '(a,es12.4)', 'factor CPU seconds: ',result%time_factorization
    print '(a,es12.4)', 'refinement CPU seconds: ',result%time_refinement
end program sparse_large_lp
