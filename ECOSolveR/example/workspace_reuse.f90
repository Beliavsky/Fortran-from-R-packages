program workspace_reuse
    use ecos
    implicit none
    integer, parameter :: n=1000
    type(ecos_dims) :: dims
    type(ecos_csc_matrix) :: g
    type(ecos_problem) :: prob
    type(ecos_workspace) :: ws
    type(ecos_result) :: first,second
    real(dp), allocatable :: c(:),h(:)
    integer :: i,ierr

    allocate(c(n),h(n),g%colptr(n+1),g%rowind(n),g%values(n))
    c=1.0_dp; h=-1.0_dp; dims%l=n
    g%nrow=n; g%ncol=n
    do i=1,n
        g%colptr(i)=i; g%rowind(i)=i; g%values(i)=-1.0_dp
    end do
    g%colptr(n+1)=n+1
    call setup_problem_csc(prob,c,g,h,dims,ierr=ierr)
    if(ierr/=0) error stop 'setup failed'
    call ecos_setup(ws,prob)
    call ecos_solve(ws,first)
    if(first%exitflag/=ECOS_OPTIMAL) error stop 'first solve failed'
    c(1)=2.0_dp
    call ecos_update(ws,c=c,ierr=ierr)
    call ecos_solve(ws,second)
    if(second%exitflag/=ECOS_OPTIMAL) error stop 'second solve failed'
    print '(a,i0)', 'first symbolic analyses: ',first%symbolic_analyses
    print '(a,i0)', 'second symbolic analyses: ',second%symbolic_analyses
    print '(a,i0)', 'second cached symbolic reuses: ',second%cached_symbolic_reuses
    print '(a,i0)', 'second cached warm starts: ',second%cached_warm_starts
end program workspace_reuse
