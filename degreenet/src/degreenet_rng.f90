! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_rng
  use degreenet_kinds, only : dp
  implicit none
  private
  public :: seed_rng, runif01, rnorm01, rexp1, rgamma_shape, rpoisson_basic

contains
  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729*i + 8191*i*i, huge(1)-1)
      if (put(i) == 0) put(i) = i
    end do
    call random_seed(put=put)
  end subroutine seed_rng

  real(dp) function runif01() result(u)
    call random_number(u)
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
    if (u >= 1.0_dp) u = 1.0_dp - epsilon(1.0_dp)
  end function runif01

  real(dp) function rnorm01() result(z)
    real(dp) :: u1, u2
    u1 = runif01(); u2 = runif01()
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
  end function rnorm01

  real(dp) function rexp1() result(x)
    x = -log(runif01())
  end function rexp1

  recursive real(dp) function rgamma_shape(shape) result(x)
    real(dp), intent(in) :: shape
    real(dp) :: d, c, z, v, u
    if (shape <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      x = rgamma_shape(shape + 1.0_dp)*runif01()**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z = rnorm01()
        v = 1.0_dp + c*z
        if (v > 0.0_dp) exit
      end do
      v = v**3
      u = runif01()
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
    end do
    x = d*v
  end function rgamma_shape

  integer function rpoisson_basic(lambda) result(k)
    real(dp), intent(in) :: lambda
    real(dp) :: l, p, z
    integer :: n
    if (lambda <= 0.0_dp) then
      k = 0
    else if (lambda < 30.0_dp) then
      l = exp(-lambda); p = 1.0_dp; n = 0
      do
        n = n + 1
        p = p*runif01()
        if (p <= l) exit
      end do
      k = n - 1
    else
      do
        z = lambda + sqrt(lambda)*rnorm01()
        if (z >= 0.0_dp) exit
      end do
      k = nint(z)
    end if
  end function rpoisson_basic
end module degreenet_rng
