! SPDX-License-Identifier: GPL-3.0-only
program test_crossdep
  use indgenerrors
  implicit none
  real(dp), parameter :: tol = 3.0e-13_dp
  real(dp) :: x(8) = [0.2_dp,-1.1_dp,0.7_dp,2.0_dp,-0.4_dp,1.3_dp,-2.2_dp,0.5_dp]
  real(dp) :: y(8) = [1.4_dp,0.1_dp,-0.8_dp,0.9_dp,2.2_dp,-1.5_dp,0.3_dp,-0.2_dp]
  real(dp) :: z(8) = [-0.5_dp,1.7_dp,0.4_dp,-1.2_dp,0.8_dp,2.1_dp,-0.9_dp,0.2_dp]
  real(dp), parameter :: es(5) = [ &
    -0.023809523809523808_dp,0.095238095238095233_dp,-0.35714285714285715_dp, &
    0.80952380952380953_dp,-0.54761904761904767_dp ]
  real(dp), parameter :: eg(5) = [ &
    -0.19807924113244763_dp,0.22404775626328199_dp,-0.29790879116899355_dp, &
    0.75147699352318587_dp,-0.55797461713249819_dp ]
  real(dp), parameter :: ee(5) = [ &
    -0.19750962141221709_dp,0.57965441978906451_dp,-0.37337429101527592_dp, &
    0.45421030490287434_dp,-0.53257713914806726_dp ]
  real(dp), parameter :: es3(9) = [ &
    0.68582765502740461_dp,0.10391328106475828_dp,-0.49878374911083972_dp, &
    -0.072739296745330792_dp,-0.50917507721731559_dp,0.24939187455541986_dp, &
    0.010391328106475828_dp,-0.29095718698132317_dp,0.25978320266189570_dp ]
  type(dependence_two_result) :: out2
  type(dependence_three_result) :: out3

  out2 = crossdep_2series(x,y,2)
  call check(maxval(abs(out2%spearman%stat-es)) < tol,'Spearman pair scores')
  call check(maxval(abs(out2%vdw%stat-eg)) < tol,'van der Waerden pair scores')
  call check(maxval(abs(out2%savage%stat-ee)) < tol,'Savage pair scores')

  out3 = crossdep_3series(x,y,z,2,1)
  call check(maxval(abs(out3%spearman%xyz%stat-es3)) < tol,'three-series Spearman scores')
  call check(out3%spearman%aggregate >= out3%spearman%xyz%aggregate,'global dependence aggregate')
  call check(out3%spearman%p_aggregate >= 0.0_dp .and. &
    out3%spearman%p_aggregate <= 1.0_dp,'global p-value range')
  print '(a)', 'test_crossdep: PASS'

contains

  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) then
      print '(a)', 'FAIL: '//message
      error stop 1
    end if
  end subroutine check

end program test_crossdep
