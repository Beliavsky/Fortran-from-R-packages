! SPDX-License-Identifier: GPL-3.0-only
program test_crosscor
  use indgenerrors
  implicit none
  real(dp), parameter :: tol = 2.0e-13_dp
  real(dp) :: x(8) = [0.2_dp,-1.1_dp,0.7_dp,2.0_dp,-0.4_dp,1.3_dp,-2.2_dp,0.5_dp]
  real(dp) :: y(8) = [1.4_dp,0.1_dp,-0.8_dp,0.9_dp,2.2_dp,-1.5_dp,0.3_dp,-0.2_dp]
  real(dp) :: z(8) = [-0.5_dp,1.7_dp,0.4_dp,-1.2_dp,0.8_dp,2.1_dp,-0.9_dp,0.2_dp]
  real(dp), parameter :: e2(5) = [ &
    -0.26344086021505381_dp,0.22759856630824379_dp,-0.22222222222222221_dp, &
    0.69534050179211471_dp,-0.56630824372759858_dp ]
  real(dp), parameter :: e3(9) = [ &
    0.79139876428316747_dp,0.00010122518793108826_dp,-0.38659923774657518_dp, &
    -0.020366507811743536_dp,-0.46148563177802571_dp,0.21101402676123543_dp, &
    -0.020346262774157318_dp,-0.36995781685069723_dp,0.19744985157846390_dp ]
  type(lag_test_result) :: out2
  type(four_lag_test_result) :: out3

  out2 = crosscor_2series(x,y,2)
  call check(out2%status == indgen_success,'crosscor_2series status')
  call check(maxval(abs(out2%stat-e2)) < tol,'crosscor_2series statistics')
  call check(abs(out2%aggregate-7.7983068048971633_dp) < 5.0e-13_dp,'crosscor_2series H')
  call check(all(out2%lags(1,:) == [-2,-1,0,1,2]),'pair lag ordering')

  out3 = crosscor_3series(x,y,z,2,1)
  call check(maxval(abs(out3%xyz%stat-e3)) < tol,'crosscor_3series triple statistics')
  call check(abs(out3%aggregate-33.991420787937656_dp) < 2.0e-12_dp,'crosscor_3series H')
  call check(all(out3%xyz%lags(:,1) == [-1,-1]),'triple lag start')
  call check(all(out3%xyz%lags(:,9) == [1,1]),'triple lag end')
  print '(a)', 'test_crosscor: PASS'

contains

  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) then
      print '(a)', 'FAIL: '//message
      error stop 1
    end if
  end subroutine check

end program test_crosscor
