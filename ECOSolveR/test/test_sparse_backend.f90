program test_sparse_backend
    use ecos
    implicit none
    integer, parameter :: nlarge=2000
    type(ecos_dims) :: d
    type(ecos_result) :: r
    type(ecos_problem) :: p
    type(ecos_csc_matrix) :: gc,ac,gl
    real(dp) :: c3(3),g3(3,3),h3(3),a2(2,3),b2(2)
    real(dp), allocatable :: cl(:),hl(:)
    integer :: i,ierr

    ! Sparse SOC: min t, x=3, y=4, (t,x,y) in Q3.
    d%l=0; d%e=0; allocate(d%q(1)); d%q=[3]
    c3=[1.0_dp,0.0_dp,0.0_dp]
    g3=0.0_dp; g3(1,1)=-1.0_dp; g3(2,2)=-1.0_dp; g3(3,3)=-1.0_dp
    h3=0.0_dp; a2=0.0_dp; a2(1,2)=1.0_dp; a2(2,3)=1.0_dp; b2=[3.0_dp,4.0_dp]
    call make_csc_matrix(g3,gc); call make_csc_matrix(a2,ac)
    call ecos_csolve(c3,gc,h3,d,r,ac,b2)
    call check(r%exitflag==ECOS_OPTIMAL,'sparse SOC status')
    call check(r%sparse_backend_used,'sparse SOC backend flag')
    call check(maxval(abs(r%x-[5.0_dp,3.0_dp,4.0_dp]))<2.0e-7_dp,'sparse SOC solution')
    call check(r%kkt_nnz>0 .and. r%ldl_nnz>0,'sparse diagnostics')

    ! Sparse exponential cone: min b subject to (1,b,1) in K_exp.
    deallocate(d%q); d%e=1
    c3=[0.0_dp,1.0_dp,0.0_dp]
    a2=0.0_dp; a2(1,1)=1.0_dp; a2(2,3)=1.0_dp; b2=[1.0_dp,1.0_dp]
    call make_csc_matrix(a2,ac)
    call ecos_csolve(c3,gc,h3,d,r,ac,b2, &
        control=ecos_control(maxit=150,feastol=1.0e-9_dp,reltol=1.0e-9_dp,abstol=1.0e-9_dp))
    call check(r%exitflag==ECOS_OPTIMAL,'sparse exp status')
    call check(abs(r%x(2)-exp(1.0_dp))<1.0e-7_dp,'sparse exp solution')

    ! Large diagonal LP.  The sparse problem object must not allocate dense G.
    d%e=0; d%l=nlarge
    allocate(cl(nlarge),hl(nlarge),gl%colptr(nlarge+1),gl%rowind(nlarge),gl%values(nlarge))
    cl=1.0_dp; hl=-1.0_dp
    gl%nrow=nlarge; gl%ncol=nlarge
    do i=1,nlarge
        gl%colptr(i)=i; gl%rowind(i)=i; gl%values(i)=-1.0_dp
    end do
    gl%colptr(nlarge+1)=nlarge+1
    call setup_problem_csc(p,cl,gl,hl,d,ierr=ierr)
    call check(ierr==0,'large sparse setup')
    call check(.not.allocated(p%gmat),'CSC setup does not densify G')
    call ecos_csolve(p,r,ecos_control(maxit=80))
    call check(r%exitflag==ECOS_OPTIMAL,'large sparse LP status')
    call check(maxval(abs(r%x-1.0_dp))<1.0e-7_dp,'large sparse LP solution')
    call check(r%kkt_nnz<=3*nlarge+10,'large sparse KKT nnz')
    call check(r%ldl_nnz<=2*nlarge,'large sparse LDL fill')
    call check(r%symbolic_analyses==1,'symbolic factorization reused')
    call check(r%numeric_factorizations>=1,'numeric refactorizations recorded')

    print '(a)', 'PASS test_sparse_backend'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(*),intent(in)::msg
        if(.not.ok) then; print '(a,1x,a)','FAIL',msg; error stop 1; end if
    end subroutine check
end program test_sparse_backend
