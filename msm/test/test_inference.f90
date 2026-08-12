program test_inference
    use msm, only : dp,delta_method
    implicit none
    real(dp) :: theta(2),cov(2,2)
    real(dp),allocatable :: est(:),ocov(:,:),se(:)
    theta=[log(2.0_dp),3.0_dp]; cov=reshape([0.04_dp,0.01_dp,0.01_dp,0.09_dp],[2,2])
    call delta_method(fun,theta,cov,est,ocov,se,1e-6_dp)
    call check(abs(est(1)-2.0_dp)<1e-14_dp,"delta estimate")
    call check(abs(ocov(1,1)-4.0_dp*0.04_dp)<1e-8_dp,"delta covariance")
    call check(abs(ocov(2,2)-(0.04_dp+2.0_dp*0.01_dp+0.09_dp))<1e-8_dp,"delta covariance linear")
    print '(a)', 'test_inference: PASS'
contains
    subroutine fun(x,y)
        real(dp),intent(in)::x(:); real(dp),allocatable,intent(out)::y(:)
        allocate(y(2)); y=[exp(x(1)),x(1)+x(2)]
    end subroutine fun
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(*),intent(in)::msg
        if(.not.ok) then; write(*,'(a)') 'FAIL: '//msg; error stop 1; end if
    end subroutine check
end program test_inference
