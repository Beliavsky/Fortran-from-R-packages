program test_lifecycle_csc
    use ecos
    implicit none
    type(ecos_dims)::d
    type(ecos_problem)::p
    type(ecos_workspace)::ws
    type(ecos_result)::r1,r2,r3
    type(ecos_csc_matrix)::gc
    real(dp)::c(2),g(2,2),h(2),a(1,2),b(1)
    real(dp), allocatable :: gd(:,:)
    integer::ierr

    c=[-1.0_dp,-1.0_dp];g=0.0_dp;g(1,1)=-1.0_dp;g(2,2)=-1.0_dp
    h=0.0_dp;a=1.0_dp;b=1.0_dp;d%l=2
    call setup_problem_dense(p,c,g,h,d,a,b,ierr=ierr)
    call check(ierr==0,'setup')
    call ecos_setup(ws,p)
    call ecos_solve(ws,r1,ierr=ierr)
    call check(ierr==0 .and. r1%exitflag==ECOS_OPTIMAL,'first solve')
    call ecos_update(ws,c=[-2.0_dp,-1.0_dp],ierr=ierr)
    call check(ierr==0,'update')
    call ecos_solve(ws,r2)
    call check(r2%exitflag==ECOS_OPTIMAL,'updated solve')
    call check(r2%x(1)>0.999999_dp,'updated optimum')
    call ecos_solve(ws,r3,control=ecos_control(maxit=1))
    call check(r3%exitflag==ECOS_MAXIT,'per-solve maxit')
    call ecos_cleanup(ws)
    call ecos_solve(ws,r3,ierr=ierr)
    call check(ierr/=0,'cleanup guard')

    call make_csc_matrix(g,gc)
    call gc%to_dense(gd)
    call check(maxval(abs(gd-g))<1.0e-15_dp,'CSC roundtrip')
    call ecos_csolve(c,gc,h,d,r1)
    call check(r1%exitflag==ECOS_DINF,'CSC unconstrained LP status')
    print '(a)', 'PASS test_lifecycle_csc'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(*),intent(in)::msg
        if(.not.ok) then; print '(a,1x,a)','FAIL',msg; error stop 1; end if
    end subroutine
end program
