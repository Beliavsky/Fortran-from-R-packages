program test_regression
    use bigstatsr
    implicit none
    type(fbm_real) :: x
    type(univ_reg_result) :: lr,gr
    real(dp), allocatable :: a(:,:),y(:),yl(:),z(:,:)
    integer :: n,i

    n=120
    allocate(a(n,3),y(n),yl(n),z(n,1))
    do i=1,n
        z(i,1)=real(i-(n+1)/2,dp)/real(n,dp)
        a(i,1)=sin(0.31_dp*real(i,dp))+0.2_dp*cos(0.07_dp*real(i,dp))
        a(i,2)=cos(0.17_dp*real(i,dp))
        a(i,3)=sin(0.11_dp*real(i,dp))
    end do
    y=1.5_dp+2.0_dp*a(:,1)+0.7_dp*z(:,1)
    yl=1.0_dp/(1.0_dp+exp(-(-0.3_dp+0.8_dp*a(:,1)+0.4_dp*z(:,1))))
    x=create_fbm('test_regression.bk',n,3)
    call fbm_from_array(x,a)

    lr=big_univ_linreg(x,y,z)
    call check(abs(lr%estim(1)-2.0_dp)<1.0e-9_dp,'linear slope')
    call check(lr%std_err(1)<1.0e-7_dp,'linear se')

    gr=big_univ_logreg(x,yl,z,tol=1.0e-10_dp,maxiter=60)
    call check(gr%converged(1),'logistic converged')
    call check(abs(gr%estim(1)-0.8_dp)<2.0e-7_dp,'logistic slope')

    call execute_command_line('rm -f test_regression.bk')
    print *, 'test_regression: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(len=*),intent(in)::msg
        if(.not.ok) then
            print *, 'FAIL: ',trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_regression
