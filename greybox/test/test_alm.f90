program test_alm
    use greybox_kinds, only: dp
    use greybox_regression, only: alm_model, alm_fit, alm_predict, point_aic, point_bic
    implicit none
    integer, parameter :: n=40
    real(dp) :: x(n,2), y(n), ylog(n), pred(n), t
    type(alm_model) :: fit
    integer :: i, failures

    failures=0
    do i=1,n
        t=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
        x(i,:)=[1.0_dp,t]
        y(i)=1.5_dp+2.25_dp*t+0.08_dp*sin(7.0_dp*t)
        ylog(i)=exp(0.4_dp+0.7_dp*t+0.03_dp*cos(5.0_dp*t))
    end do

    call alm_fit(x,y,'dnorm',fit,max_iter=800,tol=1.0e-7_dp)
    call check_true(fit%converged,'normal converged')
    call check_close(fit%beta(1),1.5_dp,0.03_dp,'normal intercept')
    call check_close(fit%beta(2),2.25_dp,0.03_dp,'normal slope')
    call alm_predict(fit,x,pred)
    call check_true(sqrt(sum((pred-y)**2)/real(n,dp))<0.08_dp,'normal prediction')
    call check_true(fit%scale>0.0_dp,'normal scale')
    call check_close(sum(point_aic(fit))/real(n,dp),fit%aic,1.0e-10_dp,'point AIC identity')
    call check_close(sum(point_bic(fit))/real(n,dp),fit%bic,1.0e-10_dp,'point BIC identity')

    call alm_fit(x,ylog,'dlnorm',fit,max_iter=900,tol=1.0e-7_dp)
    call check_true(fit%converged,'lognormal converged')
    call check_close(fit%beta(1),0.4_dp,0.04_dp,'lognormal intercept')
    call check_close(fit%beta(2),0.7_dp,0.04_dp,'lognormal slope')
    call alm_predict(fit,x,pred)
    call check_true(all(pred>0.0_dp),'lognormal positive prediction')

    call alm_fit(x,y,'dlaplace',fit,max_iter=900,tol=1.0e-7_dp)
    call check_true(fit%converged,'laplace converged')
    call check_close(fit%beta(2),2.25_dp,0.08_dp,'laplace slope')

    if(failures/=0)then
        write(*,'(a,i0)')'test_alm: FAIL ',failures
        error stop 1
    end if
    write(*,'(a)')'test_alm: PASS'
contains
    subroutine check_close(a,b,tol,name)
        real(dp),intent(in)::a,b,tol
        character(len=*),intent(in)::name
        if(abs(a-b)>tol)then
            failures=failures+1
            write(*,'(a,2f14.7)')trim(name)//': ',a,b
        end if
    end subroutine check_close
    subroutine check_true(ok,name)
        logical,intent(in)::ok
        character(len=*),intent(in)::name
        if(.not.ok)then;failures=failures+1;write(*,'(a)')trim(name)//': failed';end if
    end subroutine check_true
end program test_alm
