! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
module greeks_rng
  use, intrinsic :: iso_fortran_env, only: int64
  use greeks_kinds, only: dp, pi
  implicit none
  private
  integer(int64), save :: state = 88172645463325252_int64
  logical, save :: has_spare = .false.
  real(dp), save :: spare = 0.0_dp
  public :: seed_rng, uniform_random, normal_random, poisson_random
  public :: student_t3_random
contains
  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    state = int(max(1, seed), int64)
    state = ieor(state, shiftl(state, 21))
    state = ieor(state, shiftr(state, 35))
    state = ieor(state, shiftl(state, 4))
    if (state == 0_int64) state = 88172645463325252_int64
    has_spare = .false.
  end subroutine seed_rng

  function next_u64() result(x)
    integer(int64) :: x
    x = state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    state = x
  end function next_u64

  function uniform_random() result(value)
    real(dp) :: value
    integer(int64) :: x
    x = iand(next_u64(), int(z'7FFFFFFFFFFFFFFF', int64))
    value = (real(x, dp) + 0.5_dp)/real(huge(1_int64), dp)
    value = min(max(value, epsilon(1.0_dp)), 1.0_dp-epsilon(1.0_dp))
  end function uniform_random

  function normal_random() result(value)
    real(dp) :: value, u1, u2, radius
    if (has_spare) then
      value = spare
      has_spare = .false.
      return
    end if
    u1 = uniform_random()
    u2 = uniform_random()
    radius = sqrt(-2.0_dp*log(u1))
    value = radius*cos(2.0_dp*pi*u2)
    spare = radius*sin(2.0_dp*pi*u2)
    has_spare = .true.
  end function normal_random

  function poisson_random(lambda) result(value)
    real(dp), intent(in) :: lambda
    integer :: value, k
    real(dp) :: limit, product, z
    if (lambda <= 0.0_dp) then
      value = 0
    else if (lambda < 30.0_dp) then
      limit = exp(-lambda)
      product = 1.0_dp
      k = 0
      do
        k = k + 1
        product = product*uniform_random()
        if (product <= limit) exit
      end do
      value = k - 1
    else
      do
        z = lambda + sqrt(lambda)*normal_random()
        if (z >= 0.0_dp) exit
      end do
      value = nint(z)
    end if
  end function poisson_random

  function student_t3_random() result(value)
    real(dp) :: value, z, chi
    z = normal_random()
    chi = normal_random()**2 + normal_random()**2 + normal_random()**2
    value = z/sqrt(max(chi/3.0_dp, tiny(1.0_dp)))
  end function student_t3_random
end module greeks_rng
