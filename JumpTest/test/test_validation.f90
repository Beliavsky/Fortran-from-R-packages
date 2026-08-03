! SPDX-License-Identifier: MIT
program test_validation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use jumptest, only : dp, i8, JT_INVALID_ARGUMENT, JT_INVALID_DIMENSION, &
    JT_DEGENERATE_SAMPLE, statp_result, adjp_result, pcombine_result, &
    simulation_result, jumptestday, pcombine, ppool, sv1f, sv2f
  implicit none

  type(statp_result) :: one
  type(adjp_result) :: pooled
  type(pcombine_result) :: combined
  type(simulation_result) :: sim
  real(dp) :: pmat(2, 2)
  character(len=4) :: one_method(1)

  call jumptestday([0.0_dp, 0.0_dp, 0.0_dp], one, 'BNS')
  call check(one%status == JT_DEGENERATE_SAMPLE, 'degenerate BNS status')
  call check(ieee_is_nan(one%pvalue), 'degenerate BNS p-value')

  call jumptestday([0.01_dp, 0.02_dp, 0.03_dp], one, 'bad')
  call check(one%status == JT_INVALID_ARGUMENT, 'invalid method')

  one_method = ['BNS ']
  call pcombine(reshape([0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp], [2, 2]), &
    one_method, combined)
  call check(combined%status == JT_INVALID_DIMENSION, 'pcombine requires two methods')

  pmat = 0.5_dp
  call ppool(pmat, pooled, 'bad')
  call check(pooled%status == JT_INVALID_ARGUMENT, 'invalid pool method')

  call sv1f(10, 2, sim, correlation=1.1_dp, seed=1_i8)
  call check(sim%status == JT_INVALID_ARGUMENT, 'invalid correlation')

  call sv2f(10, 2, sim, r1=0.9_dp, r2=0.9_dp, seed=1_i8)
  call check(sim%status == JT_INVALID_ARGUMENT, 'non-positive covariance')

  call sv1f(0, 2, sim)
  call check(sim%status == JT_INVALID_DIMENSION, 'invalid simulation size')

  print '(a)', 'test_validation: PASS'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      print '(a,1x,a)', 'FAIL:', message
      error stop 1
    end if
  end subroutine check

end program test_validation
