program test_xreg_parity
  use tsa, only : dp, arimax_result, arimax_fit
  use tseries_linalg, only : solve_linear, invert_matrix
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  implicit none
  integer, parameter :: n = 80
  real(dp) :: y(n), xreg(n,2), t, beta_ref(2), xty(2), xtx(2,2)
  real(dp) :: inv(2,2), resid(n), sigma2, cov_ref(2,2), nanv, fixed(2)
  real(dp) :: ymiss(n), xmiss(n,2), xc(n,2), yc(n), beta_missing(2), xtxm(2,2), xtym(2)
  type(arimax_result) :: fit
  integer :: i, status, nc

  do i = 1, n
    t = real(i,dp)/real(n,dp)
    xreg(i,1) = t + 0.15_dp*sin(0.37_dp*real(i,dp))
    xreg(i,2) = 1.0002_dp*xreg(i,1) + 0.012_dp*cos(0.29_dp*real(i,dp))
    y(i) = 2.25_dp*xreg(i,1) - 1.70_dp*xreg(i,2) &
      + 0.04_dp*sin(0.83_dp*real(i,dp))
  end do

  xtx = matmul(transpose(xreg),xreg)
  xty = matmul(transpose(xreg),y)
  call solve_linear(xtx,xty,beta_ref,status)
  call check(status == 0,'reference OLS solve')
  resid = y-matmul(xreg,beta_ref)
  sigma2 = sum(resid*resid)/real(n,dp)
  call invert_matrix(xtx,inv,status)
  call check(status == 0,'reference covariance inverse')
  cov_ref = sigma2*inv

  fit = arimax_fit(y,0,0,0,xreg,include_mean=.false.,method='ML', &
    max_iterations=2000,tolerance=1.0e-11_dp)
  call check(fit%status == 0,'rotated xreg fit status')
  call close(fit%regression(1),beta_ref(1),3.0e-4_dp,'xreg beta1')
  call close(fit%regression(2),beta_ref(2),3.0e-4_dp,'xreg beta2')
  call close(fit%coefficients(1),beta_ref(1),3.0e-4_dp,'coef beta1')
  call close(fit%coefficients(2),beta_ref(2),3.0e-4_dp,'coef beta2')
  call close(fit%sigma2,sigma2,3.0e-4_dp,'xreg sigma2')
  call close(fit%covariance(1,1),cov_ref(1,1),2.0e-2_dp,'xreg cov11')
  call close(fit%covariance(1,2),cov_ref(1,2),2.0e-2_dp,'xreg cov12')
  call close(fit%covariance(2,2),cov_ref(2,2),2.0e-2_dp,'xreg cov22')

  nanv = ieee_value(0.0_dp,ieee_quiet_nan)
  fixed = [nanv,-1.70_dp]
  fit = arimax_fit(y,0,0,0,xreg,include_mean=.false.,method='ML',fixed=fixed, &
    max_iterations=2000,tolerance=1.0e-11_dp)
  call check(fit%status == 0,'fixed xreg fit status')
  call close(fit%regression(2),-1.70_dp,1.0e-12_dp,'fixed xreg beta2')

  ymiss = y
  xmiss = xreg
  ymiss(13) = nanv
  xmiss(27,1) = nanv
  xmiss(51,2) = nanv
  nc = 0
  do i = 1, n
    if (ieee_is_finite(ymiss(i)) .and. all(ieee_is_finite(xmiss(i,:)))) then
      nc = nc+1
      xc(nc,:) = xmiss(i,:)
      yc(nc) = ymiss(i)
    end if
  end do
  xtxm = matmul(transpose(xc(1:nc,:)),xc(1:nc,:))
  xtym = matmul(transpose(xc(1:nc,:)),yc(1:nc))
  call solve_linear(xtxm,xtym,beta_missing,status)
  call check(status == 0,'missing xreg reference solve')

  fit = arimax_fit(ymiss,0,0,0,xmiss,include_mean=.false.,method='ML', &
    max_iterations=2000,tolerance=1.0e-11_dp)
  call check(fit%status == 0,'missing xreg fit status')
  call close(fit%regression(1),beta_missing(1),5.0e-4_dp,'missing xreg beta1')
  call close(fit%regression(2),beta_missing(2),5.0e-4_dp,'missing xreg beta2')

  print '(a)', 'test_xreg_parity: PASS'
contains
  subroutine close(a,b,tol,msg)
    real(dp), intent(in) :: a,b,tol
    character(len=*), intent(in) :: msg
    if (abs(a-b) > tol*max(1.0_dp,abs(b))) then
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
end program test_xreg_parity
