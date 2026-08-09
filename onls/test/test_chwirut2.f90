program test_chwirut2
    use onls
    implicit none
    real(dp),parameter::x(54)=[ &
    0.5_dp,1.0_dp,1.75_dp,3.75_dp,5.75_dp,0.875_dp,2.25_dp,3.25_dp,5.25_dp,0.75_dp, &
    1.75_dp,2.75_dp,4.75_dp,0.625_dp,1.25_dp,2.25_dp,4.25_dp,0.5_dp,3.0_dp,0.75_dp, &
    3.0_dp,1.5_dp,6.0_dp,3.0_dp,6.0_dp,1.5_dp,3.0_dp,0.5_dp,2.0_dp,4.0_dp,0.75_dp, &
    2.0_dp,5.0_dp,0.75_dp,2.25_dp,3.75_dp,5.75_dp,3.0_dp,0.75_dp,2.5_dp,4.0_dp, &
    0.75_dp,2.5_dp,4.0_dp,0.75_dp,2.5_dp,4.0_dp,0.5_dp,6.0_dp,3.0_dp,0.5_dp,2.75_dp, &
    0.5_dp,1.75_dp]
    real(dp),parameter::y(54)=[ &
    92.9_dp,57.1_dp,31.05_dp,11.5875_dp,8.025_dp,63.6_dp,21.4_dp,14.25_dp,8.475_dp, &
    63.8_dp,26.8_dp,16.4625_dp,7.125_dp,67.3_dp,41.0_dp,21.15_dp,8.175_dp,81.5_dp, &
    13.12_dp,59.9_dp,14.62_dp,32.9_dp,5.44_dp,12.56_dp,5.44_dp,32.0_dp,13.95_dp, &
    75.8_dp,20.0_dp,10.42_dp,59.5_dp,21.67_dp,8.55_dp,62.0_dp,20.2_dp,7.76_dp,3.75_dp, &
    11.81_dp,54.7_dp,23.7_dp,11.55_dp,61.3_dp,17.7_dp,8.74_dp,59.2_dp,16.3_dp,8.62_dp, &
    81.0_dp,4.87_dp,14.62_dp,81.7_dp,17.17_dp,81.3_dp,28.9_dp]
    type(onls_result)::fit
    call fit_onls(chwirut,x,y,[0.1_dp,0.005_dp,0.01_dp],fit)
    if (.not.fit%converged) error stop 'Chwirut2 convergence'
    if (maxval(abs(fit%par_onls-[0.15406988_dp,0.0049969313_dp,0.012613780_dp])) > 2.0e-5_dp) then
        print *, fit%par_onls
        error stop 'Chwirut2 parameters'
    end if
    if (abs(fit%rss_orthogonal-5.5723405_dp) > 2.0e-4_dp) error stop 'Chwirut2 RSS'
contains
    subroutine chwirut(xx,par,yy,ierr)
        real(dp),intent(in)::xx(:),par(:)
        real(dp),intent(out)::yy(:)
        integer,intent(out)::ierr
        if(par(2)+par(3)*minval(xx)<=0.0_dp) then
            ierr=1; yy=0.0_dp; return
        end if
        yy=exp(-par(1)*xx)/(par(2)+par(3)*xx)
        ierr=0
    end subroutine chwirut
end program test_chwirut2
