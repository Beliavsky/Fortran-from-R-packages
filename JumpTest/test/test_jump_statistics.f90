! SPDX-License-Identifier: MIT
program test_jump_statistics
  use jumptest, only : dp, JT_SUCCESS, statp_result, adjp_result, pcombine_result, &
    bns_statistic, amin_statistic, amed_statistic, jumptestday, jumptestperiod, pcombine
  implicit none

  real(dp), parameter :: ret(8) = [0.01_dp, -0.02_dp, 0.015_dp, 0.04_dp, &
    -0.03_dp, 0.005_dp, -0.01_dp, 0.025_dp]
  real(dp) :: retmat(8, 3), value
  integer :: status
  type(statp_result) :: one
  type(adjp_result) :: many
  type(pcombine_result) :: combined
  character(len=4) :: methods(3)

  value = bns_statistic(ret, status)
  call check(status == JT_SUCCESS, 'BNS status')
  call check_close(value, -0.8769624565379169_dp, 2.0e-13_dp, 'BNS statistic')
  call check_close(amin_statistic(ret), -0.5591021810675295_dp, 2.0e-13_dp, &
    'Amin statistic')
  call check_close(amed_statistic(ret), -0.7209645560923156_dp, 2.0e-13_dp, &
    'Amed statistic')

  call jumptestday(ret, one, 'BNS')
  call check(one%status == JT_SUCCESS, 'day status')
  call check_close(one%pvalue, 0.809746484788029_dp, 2.0e-10_dp, 'day p-value')

  retmat(:, 1) = ret
  retmat(:, 2) = 0.5_dp*ret
  retmat(:, 3) = -2.0_dp*ret
  call jumptestperiod(retmat, many, 'Amed')
  call check(many%status == JT_SUCCESS, 'period status')
  call check(maxval(abs(many%stat - many%stat(1))) < 1.0e-12_dp, &
    'scale-invariant period statistics')
  call check(all(many%adjp >= many%pvalue), 'BH values dominate p-values')

  methods = ['BNS ', 'Amed', 'Amin']
  call pcombine(retmat, methods, combined)
  call check(combined%status == JT_SUCCESS, 'pcombine status')
  call check(all(shape(combined%pvalue) == [3, 3]), 'pcombine shape')
  call check_close(combined%pvalue(1, 1), 0.809746484788029_dp, 2.0e-10_dp, &
    'pcombine BNS')

  print '(a)', 'test_jump_statistics: PASS'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      print '(a,1x,a)', 'FAIL:', message
      error stop 1
    end if
  end subroutine check

  subroutine check_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    call check(abs(actual - expected) <= tolerance, message)
  end subroutine check_close

end program test_jump_statistics
