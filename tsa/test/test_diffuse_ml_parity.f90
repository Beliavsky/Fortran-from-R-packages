program test_diffuse_ml_parity
  use tsa, only : dp, arima_fit, arimax_result
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  type(arimax_result) :: fit
  real(dp) :: x1(6), x2(7), x3(8), x4(12), nanv
  real(dp) :: fixed11(2)
  real(dp), parameter :: tol = 5.0e-9_dp

  nanv = ieee_value(0.0_dp,ieee_quiet_nan)

  x1 = [1.0_dp,2.0_dp,4.0_dp,7.0_dp,11.0_dp,16.0_dp]
  fit = arima_fit(x1,0,1,0,include_mean=.false.,method='ML')
  call check(fit%status == 0,'rw status')
  call close(fit%sigma2,11.000000000001723_dp,tol,'rw sigma2')
  call close(fit%loglik,-13.089430848015374_dp,tol,'rw loglik')
  call close(fit%residuals(1),9.999995000004e-4_dp,tol,'rw diffuse residual')
  call close(fit%residuals(6),5.0_dp,tol,'rw final residual')
  call check(fit%n_cond == 0,'ML n.cond')

  x2 = [1.0_dp,2.0_dp,nanv,7.0_dp,11.0_dp,16.0_dp,22.0_dp]
  fit = arima_fit(x2,0,1,0,include_mean=.false.,method='ML')
  call check(fit%status == 0,'missing rw status')
  call close(fit%sigma2,18.10000000000172_dp,tol,'missing rw sigma2')
  call close(fit%loglik,-14.681046101978717_dp,tol,'missing rw loglik')
  call close(fit%residuals(4),3.535533905933_dp,tol,'missing prediction residual')

  x3 = [1.0_dp,1.5_dp,2.2_dp,3.4_dp,5.1_dp,7.5_dp,10.8_dp,15.2_dp]
  fixed11 = [0.4_dp,0.2_dp]
  fit = arima_fit(x3,1,1,1,include_mean=.false.,method='ML',fixed=fixed11)
  call check(fit%status == 0,'arma11 integrated status')
  call close(fit%sigma2,2.2481954410387095_dp,2.0e-8_dp,'arma11 integrated sigma2')
  call close(fit%loglik,-12.952565788848862_dp,2.0e-8_dp,'arma11 integrated loglik')

  x4 = [1.0_dp,2.0_dp,3.0_dp,5.0_dp,2.0_dp,4.0_dp,7.0_dp,9.0_dp, &
        3.0_dp,6.0_dp,10.0_dp,13.0_dp]
  fit = arima_fit(x4,0,0,0,include_mean=.false.,seasonal_d=1,period=4,method='ML')
  call check(fit%status == 0,'seasonal diffuse status')
  call close(fit%sigma2,8.375000000039842_dp,2.0e-8_dp,'seasonal diffuse sigma2')
  call close(fit%loglik,-19.8525125764837_dp,2.0e-8_dp,'seasonal diffuse loglik')

  print '(a)', 'test_diffuse_ml_parity: PASS'
contains
  subroutine close(a,b,t,msg)
    real(dp), intent(in) :: a,b,t
    character(len=*), intent(in) :: msg
    if (abs(a-b) > t*max(1.0_dp,abs(b))) then
      print '(a,2es24.14)', trim(msg)//' FAIL: ',a,b
      error stop 1
    end if
  end subroutine close
  subroutine check(ok,msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg
    if (.not. ok) then
      print '(a)', trim(msg)//' FAIL'
      error stop 1
    end if
  end subroutine check
end program test_diffuse_ml_parity
