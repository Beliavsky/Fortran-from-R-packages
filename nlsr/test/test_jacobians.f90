program test_jacobians
    use nlsr
    implicit none
    real(dp) :: p(2),r(2),j(2,2),truth(2,2)
    integer :: method,feval,ierr,mask(2)
    p=[0.7_dp,-0.3_dp]; mask=1
    call res(p,r,ierr)
    truth(1,:)=[2.0_dp*p(1),1.0_dp]
    truth(2,:)=[cos(p(1)),-2.0_dp]
    do method=jac_forward,jac_richardson
        feval=1
        call numerical_jacobian(res,p,r,j,method,1.0e-5_dp,mask,feval,ierr)
        if (ierr/=0) error stop 'jacobian error'
        if (maxval(abs(j-truth))>2.0e-4_dp) error stop 'jacobian mismatch'
    end do
    mask=[1,0]; feval=1
    call numerical_jacobian(res,p,r,j,jac_central,1.0e-5_dp,mask,feval,ierr)
    if (maxval(abs(j(:,2)))>1.0e-15_dp) error stop 'mask ignored'
    print *, 'PASS test_jacobians'
contains
    subroutine res(x,y,ierr)
        real(dp),intent(in)::x(:); real(dp),intent(out)::y(:); integer,intent(out)::ierr
        y(1)=x(1)*x(1)+x(2); y(2)=sin(x(1))-2.0_dp*x(2); ierr=0
    end subroutine
end program
