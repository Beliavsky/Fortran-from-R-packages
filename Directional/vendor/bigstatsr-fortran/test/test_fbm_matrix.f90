program test_fbm_matrix
    use bigstatsr
    implicit none
    type(fbm_real) :: x,xt
    type(colstats_result) :: st
    real(dp) :: a(4,3),v(3),b(3,2)
    real(dp), allocatable :: got(:),cm(:,:),k(:,:),kt(:,:),arr(:,:)

    a = reshape([ &
        1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
        2.0_dp,0.0_dp,-1.0_dp,1.0_dp, &
        5.0_dp,4.0_dp,3.0_dp,2.0_dp], [4,3])
    x=create_fbm('test_fbm_matrix.bk',4,3)
    call fbm_from_array(x,a)
    arr=x%to_array()
    call check(maxval(abs(arr-a))<1.0e-13_dp,'roundtrip')

    st=big_colstats(x)
    call check(maxval(abs(st%sum-[10.0_dp,2.0_dp,14.0_dp]))<1.0e-13_dp,'col sums')
    call check(abs(st%var(1)-5.0_dp/3.0_dp)<1.0e-13_dp,'col variance')

    v=[0.5_dp,-1.0_dp,2.0_dp]
    got=big_prod_vec(x,v)
    call check(maxval(abs(got-matmul(a,v)))<1.0e-12_dp,'prod vec')
    got=big_cprod_vec(x,[1.0_dp,-2.0_dp,0.5_dp,3.0_dp])
    call check(maxval(abs(got-matmul(transpose(a),[1.0_dp,-2.0_dp,0.5_dp,3.0_dp])))<1.0e-12_dp,'cprod vec')

    b=reshape([1.0_dp,0.0_dp,-1.0_dp,2.0_dp,0.5_dp,1.0_dp],[3,2])
    cm=big_prod_mat(x,b)
    call check(maxval(abs(cm-matmul(a,b)))<1.0e-12_dp,'prod mat')

    k=big_crossprod_self(x)
    call check(maxval(abs(k-matmul(transpose(a),a)))<1.0e-12_dp,'crossprod self')
    kt=big_tcrossprod_self(x)
    call check(maxval(abs(kt-matmul(a,transpose(a))))<1.0e-12_dp,'tcrossprod self')

    xt=create_fbm('test_fbm_transpose.bk',3,4)
    call fbm_transpose(x,xt,2)
    arr=xt%to_array()
    call check(maxval(abs(arr-transpose(a)))<1.0e-12_dp,'transpose')

    call execute_command_line('rm -f test_fbm_matrix.bk test_fbm_transpose.bk')
    print *, 'test_fbm_matrix: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(len=*),intent(in)::msg
        if(.not.ok) then
            print *, 'FAIL: ',trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_fbm_matrix
