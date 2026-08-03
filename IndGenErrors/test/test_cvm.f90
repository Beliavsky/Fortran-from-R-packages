! SPDX-License-Identifier: GPL-3.0-only
program test_cvm
  use indgenerrors
  implicit none
  real(dp), parameter :: tol = 3.0e-13_dp
  real(dp) :: x(8) = [0.2_dp,-1.1_dp,0.7_dp,2.0_dp,-0.4_dp,1.3_dp,-2.2_dp,0.5_dp]
  real(dp) :: y(8) = [1.4_dp,0.1_dp,-0.8_dp,0.9_dp,2.2_dp,-1.5_dp,0.3_dp,-0.2_dp]
  real(dp) :: z(8) = [-0.5_dp,1.7_dp,0.4_dp,-1.2_dp,0.8_dp,2.1_dp,-0.9_dp,0.2_dp]
  real(dp) :: x80(80), y80(80)
  real(dp), parameter :: ecvm2(5) = [ &
    -0.93750000000000044_dp,-1.1458333333333348_dp,-0.13888888888888826_dp, &
    2.4652777777777777_dp,0.38194444444444364_dp ]
  real(dp), parameter :: ep2(5) = [ &
    0.75719999999999998_dp,0.85739999999999994_dp,0.39300000000000002_dp, &
    0.050699999999999967_dp,0.25539999999999996_dp ]
  real(dp), parameter :: ecvm3(9) = [ &
    -0.16012690698190665_dp,-1.4502922718075424_dp,-0.86468529770228841_dp, &
    -1.9169478292976689_dp,-1.1529137302697197_dp,-1.3587911821035958_dp, &
    -2.0999500087055614_dp,-1.0751378040213646_dp,-1.6790449960674101_dp ]
  type(cvm_test_result) :: out2
  type(cvm_three_result) :: out3
  integer :: i

  out2 = cvm_2series(x,y,2)
  call check(maxval(abs(out2%cvm-ecvm2)) < tol,'cvm_2series statistics')
  call check(maxval(abs(out2%p_cvm-ep2)) < tol,'cvm_2series empirical p-values')
  call check(abs(out2%wstat-0.16319444444444439_dp) < tol,'cvm_2series W')
  call check(abs(out2%fstat-11.425355910693609_dp) < 2.0e-12_dp,'cvm_2series F')
  call check(abs(out2%p_wstat-0.21154826099527269_dp) < 2.0e-12_dp,'cvm W p-value')
  call check(abs(out2%p_fstat-0.32532407469378788_dp) < 2.0e-12_dp,'cvm F p-value')

  out3 = cvm_3series(x,y,z,2,1)
  call check(maxval(abs(out3%xyz%cvm-ecvm3)) < 5.0e-13_dp,'cvm_3series triple statistics')
  call check(abs(out3%wstat-0.89546153268702267_dp) < 3.0e-12_dp,'cvm_3series W')
  call check(abs(out3%fstat-36.951034637505437_dp) < 4.0e-12_dp,'cvm_3series F')
  call check(abs(out3%p_wstat-0.18653836762270537_dp) < 3.0e-12_dp,'cvm_3series W p-value')
  call check(abs(out3%p_fstat-0.87667598758389897_dp) < 3.0e-12_dp,'cvm_3series F p-value')

  do i = 1, 80
    x80(i) = sin(0.17_dp*real(i,dp))+0.013_dp*real(i,dp)
    y80(i) = cos(0.23_dp*real(i,dp))+0.07_dp*real(mod(i-1,4),dp)
  end do
  out2 = cvm_2series(x80,y80,1)
  call check(maxval(abs(out2%cvm-[3.6585562414265875_dp, &
    2.0269590192043929_dp,1.3310828189300423_dp])) < 2.0e-11_dp, &
    'n=80 F100 statistics')
  call check(maxval(abs(out2%p_cvm-[0.026399999999999979_dp, &
    0.078200000000000047_dp,0.12660000000000005_dp])) < tol, &
    'n=80 F100 p-values')
  call check(abs(out2%p_wstat-0.0066572451382601594_dp) < 3.0e-12_dp, &
    'n=80 W p-value')
  print '(a)', 'test_cvm: PASS'

contains

  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) then
      print '(a)', 'FAIL: '//message
      error stop 1
    end if
  end subroutine check

end program test_cvm
