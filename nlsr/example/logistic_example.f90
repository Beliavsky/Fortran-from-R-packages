program logistic_example
    use nlsr
    implicit none
    real(dp)::x(10),y(10),start(3)
    type(nlsr_result)::fit
    logical::ok
    integer::i
    do i=1,10; x(i)=real(i,dp); end do
    call logistic_value(x,80.0_dp,5.0_dp,1.4_dp,y)
    call logistic_initial(x,y,start,ok)
    if (.not.ok) error stop 'could not initialize logistic model'
    call nlfb(start,10,res,fit,jacfn=jac)
    print '(a,3f16.8)', 'coefficients:',fit%coefficients
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
