program test_all_masked
    use nlsr
    implicit none
    type(nlsr_result) :: fit
    real(dp) :: start(2),lo(2),up(2)
    start=[1.0_dp,2.0_dp]; lo=start; up=start
    call nlfb(start,2,res,fit,lower=lo,upper=up)
    if (.not.fit%converged) error stop 'all masked should return'
    if (fit%feval/=1) error stop 'all masked evaluation count'
    if (maxval(abs(fit%coefficients-start))>0.0_dp) error stop 'all masked moved'
    print *, 'PASS test_all_masked'
contains
    subroutine res(p,r,ierr)
        real(dp),intent(in)::p(:); real(dp),intent(out)::r(:); integer,intent(out)::ierr
        r=[p(1)-1.0_dp,p(2)-3.0_dp]; ierr=0
    end subroutine
end program
