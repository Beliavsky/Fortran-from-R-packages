program test_arma_forecast
    use ltsa, only : dp, ltsa_error, forecast_result, tacvf_arma, prediction_variance, trench_forecast
    implicit none
    type(ltsa_error) :: error
    type(forecast_result) :: forecast_update, forecast_direct
    real(dp), allocatable :: acvf(:), variances(:), exact_variances(:)
    real(dp) :: phi(1), theta(1), sigma2, expected
    real(dp) :: z(6)
    integer :: h, j

    phi = 0.75_dp
    sigma2 = 1.7_dp
    call tacvf_arma(phi,[real(dp)::],8,sigma2,acvf,error)
    call assert_true(error%ok(),'AR tacvf failed')
    do h=0,8
        expected = sigma2*phi(1)**h/(1.0_dp-phi(1)**2)
        call assert_close(acvf(h+1),expected,2.0e-11_dp,'AR autocovariance')
    end do

    theta = 0.4_dp
    call tacvf_arma([real(dp)::],theta,4,sigma2,acvf,error)
    call assert_true(error%ok(),'MA tacvf failed')
    call assert_close(acvf(1),sigma2*(1.0_dp+theta(1)**2),1.0e-12_dp,'MA variance')
    call assert_close(acvf(2),-sigma2*theta(1),1.0e-12_dp,'MA lag one')
    call assert_close(maxval(abs(acvf(3:))),0.0_dp,1.0e-12_dp,'MA higher lags')

    call tacvf_arma(phi,[real(dp)::],20,sigma2,acvf,error)
    call prediction_variance(acvf,5,variances,error)
    call assert_true(error%ok(),'DL prediction variance failed')
    do h=1,5
        expected = sigma2*sum([(phi(1)**(2*(j-1)),j=1,h)])
        call assert_close(variances(h),expected,2.0e-10_dp,'AR prediction variance')
    end do
    call prediction_variance(acvf,5,exact_variances,error,use_durbin_levinson=.false.)
    call assert_true(error%ok(),'Trench prediction variance failed')
    call assert_vector_close(variances,exact_variances,2.0e-9_dp,'prediction methods')

    z = [0.4_dp,0.1_dp,-0.3_dp,0.7_dp,0.5_dp,-0.2_dp]
    forecast_update = trench_forecast(z,acvf,0.0_dp,4,3,update_algorithm=.true.)
    forecast_direct = trench_forecast(z,acvf,0.0_dp,4,3,update_algorithm=.false.)
    call assert_true(forecast_update%error%ok(),'updated forecast failed')
    call assert_true(forecast_direct%error%ok(),'direct forecast failed')
    call assert_matrix_close(forecast_update%forecasts,forecast_direct%forecasts,2.0e-10_dp,'forecast methods')
    call assert_matrix_close(forecast_update%sd_forecasts,forecast_direct%sd_forecasts,2.0e-10_dp,'forecast sd methods')
    do h=1,3
        call assert_close(forecast_update%forecasts(1,h),phi(1)**h*z(4),2.0e-10_dp,'AR forecast')
    end do

    print '(a)', 'test_arma_forecast: PASS'
contains
    subroutine assert_true(condition,message)
        logical,intent(in)::condition
        character(len=*),intent(in)::message
        if(.not.condition) error stop message
    end subroutine assert_true
    subroutine assert_close(actual,expected_value,tolerance,message)
        real(dp),intent(in)::actual,expected_value,tolerance
        character(len=*),intent(in)::message
        if(abs(actual-expected_value)>tolerance) then
            print *,message,actual,expected_value
            error stop
        end if
    end subroutine assert_close
    subroutine assert_vector_close(actual,expected_value,tolerance,message)
        real(dp),intent(in)::actual(:),expected_value(:),tolerance
        character(len=*),intent(in)::message
        if(maxval(abs(actual-expected_value))>tolerance) then
            print *,message,maxval(abs(actual-expected_value))
            error stop
        end if
    end subroutine assert_vector_close
    subroutine assert_matrix_close(actual,expected_value,tolerance,message)
        real(dp),intent(in)::actual(:,:),expected_value(:,:),tolerance
        character(len=*),intent(in)::message
        if(maxval(abs(actual-expected_value))>tolerance) then
            print *,message,maxval(abs(actual-expected_value))
            error stop
        end if
    end subroutine assert_matrix_close
end program test_arma_forecast
