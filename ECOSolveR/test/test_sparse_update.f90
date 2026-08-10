program test_sparse_update
    use ecos
    implicit none
    type(ecos_dims) :: dims
    type(ecos_csc_matrix) :: g,g2
    type(ecos_problem) :: prob
    type(ecos_workspace) :: ws
    type(ecos_result) :: r
    real(dp) :: c(2),h(2)
    integer :: ierr

    c=1.0_dp; h=-1.0_dp; dims%l=2
    g%nrow=2; g%ncol=2
    allocate(g%colptr(3),g%rowind(2),g%values(2))
    g%colptr=[1,2,3]; g%rowind=[1,2]; g%values=[-1.0_dp,-1.0_dp]
    call setup_problem_csc(prob,c,g,h,dims,ierr=ierr)
    call check(ierr==0,'sparse setup')
    call ecos_setup(ws,prob)
    call ecos_solve(ws,r)
    call check(r%exitflag==ECOS_OPTIMAL .and. maxval(abs(r%x-1.0_dp))<1.0e-7_dp,'first solve')

    call ecos_update(ws,h=[-2.0_dp,-2.0_dp],ierr=ierr)
    call check(ierr==0,'sparse h update')
    call ecos_solve(ws,r)
    call check(r%exitflag==ECOS_OPTIMAL .and. maxval(abs(r%x-2.0_dp))<1.0e-7_dp,'updated h solve')

    g2=g; g2%values=[-2.0_dp,-2.0_dp]
    call ecos_update(ws,g_csc=g2,h=[-2.0_dp,-2.0_dp],ierr=ierr)
    call check(ierr==0,'sparse matrix update')
    call ecos_solve(ws,r)
    call check(r%exitflag==ECOS_OPTIMAL .and. maxval(abs(r%x-1.0_dp))<1.0e-7_dp,'updated G solve')
    call check(r%sparse_backend_used,'updated solve stays sparse')

    call ecos_cleanup(ws)
    print '(a)', 'PASS test_sparse_update'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(*),intent(in)::msg
        if(.not.ok) then; print '(a,1x,a)','FAIL',msg; error stop 1; end if
    end subroutine check
end program test_sparse_update
