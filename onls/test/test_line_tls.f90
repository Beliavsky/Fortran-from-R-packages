program test_line_tls
    use onls
    implicit none
    real(dp), parameter :: x(6)=[-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp,3.0_dp]
    real(dp), parameter :: y(6)=[-4.2_dp,-0.7_dp,2.3_dp,5.4_dp,7.7_dp,11.2_dp]
    real(dp) :: start(2)
    type(onls_result) :: fit
    start=[0.0_dp,1.0_dp]
    call fit_onls(line_model,x,y,start,fit)
    if (.not.fit%converged) error stop 'line TLS did not converge'
    if (maxval(abs(fit%par_onls-[2.108728014371171_dp,3.015877304590991_dp]))>5.0e-5_dp) then
        print *, fit%par_onls
        error stop 'line TLS parameters'
    end if
    if (.not.all(fit%ortho)) error stop 'orthogonality check'
contains
    subroutine line_model(xx,par,yy,ierr)
        real(dp),intent(in)::xx(:),par(:)
        real(dp),intent(out)::yy(:)
        integer,intent(out)::ierr
        yy=par(1)+par(2)*xx
        ierr=0
    end subroutine line_model
end program test_line_tls
