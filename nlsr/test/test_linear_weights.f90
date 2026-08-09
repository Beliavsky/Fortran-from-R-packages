program test_linear_weights
    use nlsr
    implicit none
    real(dp),parameter :: t(6)=[1._dp,2._dp,3._dp,4._dp,5._dp,6._dp]
    real(dp),parameter :: y(6)=[2.0_dp,4.1_dp,5.9_dp,8.2_dp,9.8_dp,12.1_dp]
    real(dp)::start(1),w(6),expected,cov(1,1),sigma
    type(nlsr_result)::fit
    logical::ok
    start=[1.0_dp]; w=[1._dp,2._dp,1._dp,3._dp,1._dp,2._dp]
    expected=sum(w*t*y)/sum(w*t*t)
    call nlfb(start,6,res,fit,weights=w)
    if (.not.fit%converged) error stop 'weighted linear did not converge'
    if (abs(fit%coefficients(1)-expected)>1.0e-8_dp) error stop 'weighted slope'
    call covariance_from_jacobian(fit%jacobian,fit%ssquares,6,1,cov,sigma,ok)
    if (.not.ok .or. sigma<=0.0_dp .or. cov(1,1)<=0.0_dp) error stop 'covariance'
    print *, 'PASS test_linear_weights',fit%coefficients(1)
contains
    subroutine res(p,r,ierr)
        real(dp),intent(in)::p(:); real(dp),intent(out)::r(:); integer,intent(out)::ierr
        r=y-p(1)*t; ierr=0
    end subroutine
end program
