! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
program test_long_memory
  use waveslim
  use waveslim_test_support
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  real(dp) :: acvs(128)
  real(dp), allocatable :: x(:), simulated(:), variances(:)
  logical, allocatable :: basis(:)
  type(estimate_result) :: fit
  integer :: i, status
  real(dp) :: hyper

  call assert_close_scalar(fdp_sdf(0.25_dp,0.0_dp), 1.0_dp, &
    1.0e-14_dp, 'FDP white-noise spectrum')
  call assert_true(spp_sdf(0.10_dp,0.20_dp,0.18_dp) > 0.0_dp, &
    'SPP spectrum positive')
  hyper = hypergeometric_2f1(1.0_dp,1.0_dp,2.0_dp,0.25_dp,status)
  call assert_true(status == 0, 'hypergeometric convergence')
  call assert_close_scalar(hyper, -log(0.75_dp)/0.25_dp, &
    1.0e-11_dp, 'hypergeometric identity')

  do i = 1, size(acvs)
    acvs(i) = 0.8_dp**real(i-1,dp)
  end do
  x = hosking_sim(128, acvs, 12345_i8)
  call assert_true(all(ieee_is_finite(x)), 'Hosking simulation finite')
  call assert_true(sum(x*x) > 0.0_dp, 'Hosking simulation nondegenerate')

  fit = fdp_mle(x, 'haar', 4)
  call assert_true(allocated(fit%estimate), 'FDP fit estimates allocated')
  call assert_true(all(ieee_is_finite(fit%estimate)), 'FDP fit finite')

  basis = find_adaptive_basis('la8', 4, 0.18_dp, 0.05_dp)
  call assert_true(size(basis) == 30, 'adaptive-basis size')
  call assert_true(count(basis) > 0, 'adaptive basis nonempty')
  variances = bandpass_var_spp(0.2_dp,0.18_dp,4,basis)
  call assert_true(all(variances >= 0.0_dp), 'bandpass variances nonnegative')

  simulated = dwpt_sim(128,'haar',0.15_dp,0.20_dp,2,.true.,0.05_dp,77_i8)
  call assert_true(size(simulated) == 128, 'DWPT simulation size')
  call assert_true(all(ieee_is_finite(simulated)), 'DWPT simulation finite')

  write(*,'(a)') 'test_long_memory: PASS'
end program test_long_memory
