program example_alm
    use greybox_kinds, only: dp
    use greybox_regression, only: alm_model, alm_fit
    implicit none
    integer, parameter :: n=20
    real(dp) :: x(n,2), y(n), t
    type(alm_model) :: model
    integer :: i
    do i=1,n
        t=real(i-1,dp)/real(n-1,dp)
        x(i,:)=[1.0_dp,t]
        y(i)=1.0_dp+2.0_dp*t+0.05_dp*sin(10.0_dp*t)
    end do
    call alm_fit(x,y,'dnorm',model)
    write(*,'(a,2f12.6)') 'coefficients: ',model%beta
    write(*,'(a,f12.6)') 'scale: ',model%scale
    write(*,'(a,f12.6)') 'AICc: ',model%aicc
end program example_alm
