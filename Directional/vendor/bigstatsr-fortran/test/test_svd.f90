program test_svd
    use bigstatsr
    implicit none
    type(fbm_real), target :: x
    type(big_svd_result) :: s1,s2
    real(dp) :: a(6,3)
    real(dp), allocatable :: rec(:,:)
    a=0.0_dp
    a(:,1)=[1.0_dp,0.0_dp,-1.0_dp,2.0_dp,-2.0_dp,0.0_dp]
    a(:,2)=[0.0_dp,1.0_dp,1.0_dp,0.0_dp,-1.0_dp,-1.0_dp]
    a(:,3)=a(:,1)+2.0_dp*a(:,2)
    x=create_fbm('test_svd.bk',6,3)
    call fbm_from_array(x,a)
    s1=big_svd(x,2)
    call check(s1%info==0,'big_svd info')
    call check(size(s1%d)==2,'big_svd rank')
    rec=matmul(s1%u,matmul(diag2(s1%d),transpose(s1%v)))
    call check(maxval(abs(rec-a))<1.0e-9_dp,'big_svd reconstruction')
    s2=big_random_svd(x,2,tol=1.0e-10_dp,maxiter=1000,ncv=3)
    call check(s2%info==0,'big_random_svd info')
    call check(maxval(abs(s2%d-s1%d))<1.0e-7_dp,'random svd values')
    call execute_command_line('rm -f test_svd.bk')
    print *, 'test_svd: PASS'
contains
    function diag2(d) result(m)
        real(dp),intent(in)::d(:)
        real(dp)::m(size(d),size(d))
        integer::j
        m=0.0_dp
        do j=1,size(d)
            m(j,j)=d(j)
        end do
    end function diag2
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(len=*),intent(in)::msg
        if(.not.ok) then
            print *, 'FAIL: ',trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_svd
