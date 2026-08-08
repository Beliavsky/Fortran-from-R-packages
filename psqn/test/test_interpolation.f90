! SPDX-License-Identifier: Apache-2.0
program test_interpolation
  use psqn_types, only : dp
  use psqn_interpolation, only : intrapolate_type
  implicit none

  type(intrapolate_type) :: inter
  real(dp) :: val
  real(dp), parameter :: truth = 0.467251416997127_dp

  call inter%init(0.0_dp, -1.0_dp, 2.5_dp, 3.75_dp)
  val = inter%get_value(-2.0_dp, 3.0_dp)
  if (abs(val - 0.5_dp) > 1.0e-8_dp) error stop "quadratic interpolation mismatch"

  call inter%init(0.0_dp, -1.0_dp, 2.5_dp, 5.3125_dp)
  call inter%update(0.4_dp, -0.2336_dp)
  val = inter%get_value(0.4_dp, 2.5_dp)
  if (abs((val - truth) / truth) > 1.0e-8_dp) error stop "cubic interpolation mismatch"

  print *, "test_interpolation: PASS"
end program test_interpolation
