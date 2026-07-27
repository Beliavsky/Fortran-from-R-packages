! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
program test_analysis
  use pbo, only : dp, pbo_result, dominance_result, selection_result, compute_pbo, &
    column_mean, dominance_curve, selection_frequencies
  implicit none
  type(pbo_result) :: result
  type(dominance_result) :: curve
  type(selection_result) :: freq
  real(dp) :: x(12,4)
  real(dp), parameter :: cdf_selected(4) = [1.0_dp/6.0_dp, 2.0_dp/6.0_dp, &
    2.0_dp/6.0_dp, 5.0_dp/6.0_dp]
  real(dp), parameter :: cdf_all(4) = [2.0_dp/24.0_dp, 5.0_dp/24.0_dp, &
    12.0_dp/24.0_dp, 18.0_dp/24.0_dp]

  call fill_data(x)
  call compute_pbo(x,4,column_mean,result)
  call assert_true(result%success,result%message)
  call selection_frequencies(result,freq)
  call assert_true(all(freq%strategy == [1,2,4,3]),'selection ordering')
  call assert_true(all(freq%frequency == [4,1,1,0]),'selection frequencies')
  call dominance_curve(result,0.25_dp,curve)
  call assert_true(curve%success,curve%message)
  call assert_true(size(curve%performance)==4,'dominance grid size')
  call assert_vector_close(curve%cdf_selected,cdf_selected,1.0e-14_dp,'selected CDF')
  call assert_vector_close(curve%cdf_all,cdf_all,1.0e-14_dp,'all-strategy CDF')
  call assert_close(curve%integrated_difference(4),-1.0_dp/96.0_dp,1.0e-14_dp, &
    'integrated difference')
  print '(a)', 'test_analysis: PASS'
contains
  subroutine fill_data(a)
    real(dp), intent(out) :: a(12,4)
    a(1,:) = [-1.2238250364546313_dp, 1.8637284581291103_dp, -1.3706617379590857_dp, 0.0_dp]
    a(2,:) = [0.12465669298947904_dp, -0.14088465208560907_dp, -1.7768836108738526_dp, 0.8414709848078965_dp]
    a(3,:) = [0.5610581130548951_dp, -1.3528630630121898_dp, 2.0292278361970335_dp, 0.9092974268256817_dp]
    a(4,:) = [-0.5593871804245065_dp, 1.5021982742122517_dp, -0.6942259005932776_dp, 0.1411200080598672_dp]
    a(5,:) = [0.9888443445192008_dp, -0.6566681331396765_dp, 0.43949387803229234_dp, -0.7568024953079282_dp]
    a(6,:) = [1.5222980607327856_dp, 0.3003014847008945_dp, 0.8574647959705144_dp, -0.9589242746631385_dp]
    a(7,:) = [0.04181073932312873_dp, 0.04948393210667501_dp, -1.2981465270318495_dp, -0.27941549819892586_dp]
    a(8,:) = [1.9247399323163303_dp, 2.21815942636784_dp, 0.9137249801744044_dp, 0.6569865987187891_dp]
    a(9,:) = [-0.7589883130180108_dp, -1.6093882869743164_dp, -1.185019286201391_dp, 0.9893582466233818_dp]
    a(10,:) = [0.9519393955576869_dp, -1.0587603195674529_dp, -0.9104931674634585_dp, 0.4121184852417566_dp]
    a(11,:) = [0.5129029184347005_dp, -0.5308116902234399_dp, 1.679074029558842_dp, -0.5440211108893698_dp]
    a(12,:) = [0.13384911099833763_dp, -1.5082144670930706_dp, 0.6359568505505235_dp, -0.9999902065507035_dp]
  end subroutine fill_data
  subroutine assert_vector_close(actual,expected,tolerance,label)
    real(dp), intent(in) :: actual(:),expected(:),tolerance
    character(len=*), intent(in) :: label
    if (size(actual)/=size(expected) .or. maxval(abs(actual-expected))>tolerance) then
      print '(a)', 'FAILED: '//trim(label)
      error stop 1
    end if
  end subroutine assert_vector_close
  subroutine assert_close(actual,expected,tolerance,label)
    real(dp), intent(in) :: actual,expected,tolerance
    character(len=*), intent(in) :: label
    if (abs(actual-expected)>tolerance) then
      print '(a,2es24.15)', 'FAILED: '//trim(label)//' actual/expected ',actual,expected
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', 'FAILED: '//trim(label)
      error stop 1
    end if
  end subroutine assert_true
end program test_analysis
