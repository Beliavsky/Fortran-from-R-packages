program test_bounds_weights
    use onls
    implicit none
    real(dp),parameter::x(6)=[-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp,3.0_dp]
    real(dp),parameter::y(6)=[-4.2_dp,-0.7_dp,2.3_dp,5.4_dp,7.7_dp,11.2_dp]
    real(dp)::start(2),lo(2),hi(2),w(6)
    type(onls_result)::fit
    type(onls_control)::ctl
    start=[0.0_dp,1.0_dp]
    lo=[-10.0_dp,0.0_dp]; hi=[10.0_dp,2.5_dp]
    w=[1.0_dp,2.0_dp,1.0_dp,2.0_dp,1.0_dp,1.0_dp]
    ctl%mimic_r_unsorted_weights=.false.
    call fit_onls(line_model,x,y,start,fit,control=ctl,lower=lo,upper=hi,weights=w)
    if (.not.fit%converged) error stop 'bounded weighted fit'
    if(abs(fit%par_onls(2)-2.5_dp)>2.0e-5_dp) error stop 'upper bound not active'
    if(.not.(vertical_loglik(fit)<=huge(1.0_dp))) error stop 'vertical loglik'
    if(.not.(orthogonal_loglik(fit)<=huge(1.0_dp))) error stop 'orthogonal loglik'
contains
    subroutine line_model(xx,par,yy,ierr)
        real(dp),intent(in)::xx(:),par(:)
        real(dp),intent(out)::yy(:)
        integer,intent(out)::ierr
        yy=par(1)+par(2)*xx
        ierr=0
    end subroutine line_model
end program test_bounds_weights
