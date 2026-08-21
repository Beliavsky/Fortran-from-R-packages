! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_rng
  use dirichletreg_kinds, only : dp
  implicit none
  private
  public :: seed_rng, random_normal, random_gamma

contains

  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)

    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729*i + 12345*i*i, huge(1)-1)
      if (put(i) == 0) put(i) = i
    end do
    call random_seed(put=put)
  end subroutine seed_rng


  real(dp) function random_normal() result(z)
    real(dp) :: u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
  end function random_normal


  recursive real(dp) function random_gamma(shape) result(x)
    real(dp), intent(in) :: shape
    real(dp) :: d, c, z, v, u

    if (shape <= 0.0_dp) then
      x = 0.0_dp
      return
    end if

    if (shape < 1.0_dp) then
      call random_number(u)
      u = max(u, tiny(1.0_dp))
      x = random_gamma(shape + 1.0_dp)*u**(1.0_dp/shape)
      return
    end if

    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z = random_normal()
        v = 1.0_dp + c*z
        if (v > 0.0_dp) exit
      end do
      v = v*v*v
      call random_number(u)
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(max(u,tiny(1.0_dp))) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
    end do
    x = d*v
  end function random_gamma

end module dirichletreg_rng
