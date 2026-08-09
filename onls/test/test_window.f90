program test_window
    use onls
    implicit none
    integer,parameter::n=40
    real(dp)::x(n),y(n),start(2)
    type(onls_result)::fit
    type(onls_control)::ctl
    integer::i
    do i=1,n
        x(i)=-2.0_dp+4.0_dp*real(i-1,dp)/real(n-1,dp)
        y(i)=1.5_dp+0.8_dp*x(i)+0.02_dp*sin(real(i,dp))
    end do
    start=[0.0_dp,0.0_dp]
    ctl%window=6
    call fit_onls(line_model,x,y,start,fit,control=ctl)
    if (.not.fit%converged) error stop 'window fit did not converge'
    if(abs(fit%par_onls(1)-1.5_dp)>0.02_dp) error stop 'window intercept'
    if(abs(fit%par_onls(2)-0.8_dp)>0.02_dp) error stop 'window slope'
contains
    subroutine line_model(xx,par,yy,ierr)
        real(dp),intent(in)::xx(:),par(:)
        real(dp),intent(out)::yy(:)
        integer,intent(out)::ierr
        yy=par(1)+par(2)*xx
        ierr=0
    end subroutine line_model
end program test_window
