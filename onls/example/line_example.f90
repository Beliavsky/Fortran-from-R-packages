program line_example
    use onls
    implicit none
    real(dp),parameter::x(6)=[-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp,3.0_dp]
    real(dp),parameter::y(6)=[-4.2_dp,-0.7_dp,2.3_dp,5.4_dp,7.7_dp,11.2_dp]
    type(onls_result)::fit
    call fit_onls(line_model,x,y,[0.0_dp,1.0_dp],fit)
    print '(a,2f14.8)','ordinary NLS: ',fit%par_nls
    print '(a,2f14.8)','orthogonal NLS: ',fit%par_onls
    print '(a,es14.5)','orthogonal RSS: ',fit%rss_orthogonal
contains
    subroutine line_model(xx,par,yy,ierr)
        real(dp),intent(in)::xx(:),par(:)
        real(dp),intent(out)::yy(:)
        integer,intent(out)::ierr
        yy=par(1)+par(2)*xx
        ierr=0
    end subroutine line_model
end program line_example
