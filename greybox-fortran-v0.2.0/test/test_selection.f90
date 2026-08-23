program test_selection
    use greybox_kinds, only: dp
    use greybox_regression, only: alm_model
    use greybox_selection, only: calm_model, stepwise_fit, calm_fit, calm_predict, recursive_lm
    implicit none
    integer, parameter :: n=36
    real(dp) :: x(n,3), y(n), pred(n), beta_time(3,n), fitted(n), t
    logical :: selected(3)
    type(alm_model) :: step_model
    type(calm_model) :: avg_model
    integer :: i, failures

    failures=0
    do i=1,n
        t=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
        x(i,1)=1.0_dp
        x(i,2)=t
        x(i,3)=sin(11.0_dp*t)
        y(i)=2.0_dp+3.0_dp*t+0.06_dp*cos(7.0_dp*t)
    end do
    call stepwise_fit(x,y,'dnorm',step_model,selected,criterion='BIC',max_iter=500,tol=1.0e-6_dp)
    call check_true(selected(1),'stepwise intercept')
    call check_true(selected(2),'stepwise signal')
    call check_true(.not.selected(3),'stepwise noise rejection')

    call calm_fit(x,y,'dnorm',avg_model,criterion='BIC',max_iter=400,tol=1.0e-6_dp)
    call check_true(avg_model%converged,'calm converged')
    call check_true(avg_model%inclusion_probability(2)>0.95_dp,'calm signal inclusion')
    call check_true(avg_model%inclusion_probability(2)>avg_model%inclusion_probability(3),'calm ranking')
    call calm_predict(avg_model,x,pred)
    call check_true(sqrt(sum((pred-y)**2)/real(n,dp))<0.20_dp,'calm prediction')

    call recursive_lm(x,y,1.0_dp,beta_time,fitted)
    call check_true(abs(beta_time(2,n)-3.0_dp)<0.1_dp,'recursive slope')

    if(failures/=0)then
        write(*,'(a,i0)')'test_selection: FAIL ',failures
        error stop 1
    end if
    write(*,'(a)')'test_selection: PASS'
contains
    subroutine check_true(ok,name)
        logical,intent(in)::ok
        character(len=*),intent(in)::name
        if(.not.ok)then;failures=failures+1;write(*,'(a)')trim(name)//': failed';end if
    end subroutine check_true
end program test_selection
