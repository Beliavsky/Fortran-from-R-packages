program test_logistic
    use nlsr
    implicit none
    real(dp) :: x(10), y(10), start(3), v(10), jm(10,3)
    type(nlsr_result) :: fit
    type(nlsr_control) :: ctl
    logical :: ok
    integer :: i
    do i=1,10; x(i)=real(i,dp); end do
    call logistic_value(x,80.0_dp,5.0_dp,1.4_dp,y)
    call logistic_initial(x,y,start,ok)
    if (.not.ok) error stop 'logistic initial failed'
    ctl=nlsr_control(); ctl%stepredn=0.5_dp
    call nlfb(start,10,res,fit,ctl,jac)
    if (.not.fit%converged) error stop 'logistic fit did not converge'
    if (maxval(abs(fit%coefficients-[80.0_dp,5.0_dp,1.4_dp]))>2.0e-5_dp) error stop 'logistic fit'
    call logistic_value(x,fit%coefficients(1),fit%coefficients(2),fit%coefficients(3),v)
    call logistic_jacobian(x,fit%coefficients(1),fit%coefficients(2),fit%coefficients(3),jm)
    if (maxval(abs(v-y))>1.0e-5_dp) error stop 'logistic value'
    print *, 'PASS test_logistic',fit%coefficients
contains
    subroutine res(p,r,ierr)
        real(dp),intent(in)::p(:); real(dp),intent(out)::r(:); integer,intent(out)::ierr
        real(dp)::yh(10); call logistic_value(x,p(1),p(2),p(3),yh); r=y-yh; ierr=0
    end subroutine
    subroutine jac(p,j,ierr)
        real(dp),intent(in)::p(:); real(dp),intent(out)::j(:,:); integer,intent(out)::ierr
        call logistic_jacobian(x,p(1),p(2),p(3),j); j=-j; ierr=0
    end subroutine
end program
