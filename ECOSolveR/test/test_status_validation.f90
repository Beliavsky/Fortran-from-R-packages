program test_status_validation
    use ecos
    implicit none
    type(ecos_dims)::d,d0
    type(ecos_result)::r
    type(ecos_problem)::p
    real(dp)::c(1),g(2,1),h(2),g0(0,1),h0(0)
    integer::ierr
    c=0.0_dp;g(:,1)=[-1.0_dp,1.0_dp];h=[-1.0_dp,0.0_dp];d%l=2
    call ecos_csolve(c,g,h,d,r,control=ecos_control(maxit=80))
    call check(r%exitflag==ECOS_PINF,'primal infeasible')
    c=-1.0_dp
    call ecos_csolve(c,g0,h0,d0,r)
    call check(r%exitflag==ECOS_DINF,'unbounded empty problem')
    d%l=1
    call setup_problem_dense(p,c,g,h,d,ierr=ierr)
    call check(ierr/=0,'dimension validation')
    print '(a)', 'PASS test_status_validation'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(*),intent(in)::msg
        if(.not.ok) then; print '(a,1x,a)','FAIL',msg; error stop 1; end if
    end subroutine
end program
