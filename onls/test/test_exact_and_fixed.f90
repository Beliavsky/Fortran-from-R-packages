program test_exact_and_fixed
    use onls
    implicit none
    integer,parameter::n=9
    real(dp)::x(n),y(n),start(2)
    logical::fixed(2)
    type(onls_result)::fit
    integer::i
    do i=1,n
        x(i)=real(i-5,dp)/2.0_dp
        y(i)=2.0_dp+3.0_dp*x(i)
    end do
    start=[0.5_dp,0.5_dp]
    call fit_onls(line_model,x,y,start,fit)
    if (.not.fit%converged) error stop 'exact line convergence'
    if(maxval(abs(fit%par_onls-[2.0_dp,3.0_dp]))>2.0e-5_dp) error stop 'exact line result'
    fixed=[.true.,.false.]
    start=[2.0_dp,1.0_dp]
    call fit_onls(line_model,x,y,start,fit,fixed=fixed)
    if(abs(fit%par_onls(1)-2.0_dp)>1.0e-14_dp) error stop 'fixed intercept changed'
    if(abs(fit%par_onls(2)-3.0_dp)>2.0e-5_dp) error stop 'fixed slope result'
contains
    subroutine line_model(xx,par,yy,ierr)
        real(dp),intent(in)::xx(:),par(:)
        real(dp),intent(out)::yy(:)
        integer,intent(out)::ierr
        yy=par(1)+par(2)*xx
        ierr=0
    end subroutine line_model
end program test_exact_and_fixed
