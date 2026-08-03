! SPDX-License-Identifier: MIT
program test_pooling
  use jumptest, only : dp, JT_SUCCESS, adjp_result, ppool
  implicit none

  real(dp), parameter :: pmat(6, 3) = reshape([ &
    0.01_dp, 0.20_dp, 0.50_dp, 0.80_dp, 0.05_dp, 0.30_dp, &
    0.02_dp, 0.15_dp, 0.40_dp, 0.70_dp, 0.10_dp, 0.35_dp, &
    0.03_dp, 0.25_dp, 0.60_dp, 0.90_dp, 0.08_dp, 0.25_dp], [6, 3])
  type(adjp_result) :: result

  call ppool(pmat, result, 'SI')
  call check(result%status == JT_SUCCESS, 'SI status')
  call check_close(result%stat(1), 3.61472675_dp, 2.0e-8_dp, 'SI statistic')
  call check_close(result%pvalue(5), 0.00619601936_dp, 2.0e-9_dp, 'SI p-value')

  call ppool(pmat, result, 'FI')
  call check_close(result%stat(1), 24.04750218_dp, 2.0e-8_dp, 'FI statistic')
  call check_close(result%pvalue(1), 0.000511854277_dp, 2.0e-10_dp, 'FI p-value')

  call ppool(pmat, result, 'MI')
  call check_close(result%pvalue(1), 0.029701_dp, 2.0e-12_dp, 'minimum p-value')
  call check_close(result%adjp(1), 0.178206_dp, 2.0e-12_dp, 'minimum BH')

  call ppool(pmat, result, 'MA')
  call check_close(result%pvalue(4), 0.729_dp, 2.0e-12_dp, 'maximum p-value')

  call ppool(pmat, result, 'SD')
  call check_close(result%stat(1), 2.11161634_dp, 3.0e-8_dp, 'dependent Stouffer')
  call check_close(result%pvalue(1), 0.01735968_dp, 3.0e-8_dp, 'dependent Stouffer p')

  call ppool(pmat, result, 'FD')
  call check_close(result%stat(1), 8.35267319_dp, 3.0e-8_dp, 'dependent Fisher')
  call check_close(result%pvalue(1), 0.01682227_dp, 3.0e-8_dp, 'dependent Fisher p')

  print '(a)', 'test_pooling: PASS'

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

end program test_pooling
