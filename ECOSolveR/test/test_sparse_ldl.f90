program test_sparse_ldl
    use ecos_types, only : dp, ecos_csc_matrix
    use ecos_sparse, only : sparse_triplet_builder, sparse_ldl_factor, triplet_to_csc, symmetric_csc_matvec
    implicit none
    integer, parameter :: n=8
    type(sparse_triplet_builder) :: tb
    type(ecos_csc_matrix) :: a
    type(sparse_ldl_factor) :: fac
    real(dp) :: b(n),x(n),ax(n),sgn(n)
    integer :: i,info,nref

    ! Symmetric quasi-definite tridiagonal matrix with alternating signs.
    call tb%init(n,n,3*n)
    do i=1,n
        if(mod(i,2)==1) then
            call tb%add(i,i,3.0_dp)
        else
            call tb%add(i,i,-3.0_dp)
        end if
        if(i<n) call tb%add(i,i+1,0.2_dp)
    end do
    call triplet_to_csc(tb,a,.true.)
    do i=1,n
        b(i)=real(i,dp)/7.0_dp
        if(mod(i,2)==1) then; sgn(i)=1.0_dp; else; sgn(i)=-1.0_dp; end if
    end do
    call fac%analyze(a,.true.,info); call check(info==0,'LDL symbolic')
    call fac%factorize(a,sgn,1.0e-10_dp,1.0e-14_dp,info); call check(info==0,'LDL numeric')
    call fac%solve_refined(a,b,x,4,1.0e-13_dp,nref,info); call check(info==0,'LDL solve')
    call symmetric_csc_matvec(a,x,ax)
    call check(maxval(abs(ax-b))<1.0e-11_dp,'LDL residual')
    call check(fac%symbolic_nnz<n*n/2,'LDL remains sparse')
    print '(a)', 'PASS test_sparse_ldl'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(*),intent(in)::msg
        if(.not.ok) then; print '(a,1x,a)','FAIL',msg; error stop 1; end if
    end subroutine check
end program test_sparse_ldl
