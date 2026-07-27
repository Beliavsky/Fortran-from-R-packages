! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2011-2025 Genaro Sucarrat
! Copyright (C) 2026 contributors to the Modern Fortran translation
!
! This file is part of betategarch-modern-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License version 2 only.

module betategarch_rng
  use betategarch_kinds, only : dp
  use betategarch_math, only : pi
  implicit none
  private

  public :: set_random_seed, random_normal, random_gamma, random_chisq, random_student_t

  real(dp), save :: normal_spare = 0.0_dp
  logical, save :: normal_has_spare = .false.

contains

  subroutine set_random_seed(seed)
    integer, intent(in) :: seed

    integer, allocatable :: put(:)
    integer :: n, i

    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729*i + 8191*i*i, huge(1) - 1)
      if (put(i) == 0) put(i) = i
    end do
    call random_seed(put=put)
    normal_spare = 0.0_dp
    normal_has_spare = .false.
  end subroutine set_random_seed

  function random_normal() result(x)
    real(dp) :: x

    real(dp) :: u1, u2, radius

    if (normal_has_spare) then
      x = normal_spare
      normal_has_spare = .false.
      return
    end if

    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    radius = sqrt(-2.0_dp*log(u1))
    x = radius*cos(2.0_dp*pi*u2)
    normal_spare = radius*sin(2.0_dp*pi*u2)
    normal_has_spare = .true.
  end function random_normal

  recursive function random_gamma(shape) result(x)
    real(dp), intent(in) :: shape
    real(dp) :: x

    real(dp) :: d, c, z, u, v

    if (shape <= 0.0_dp) then
      x = 0.0_dp
      return
    end if

    if (shape < 1.0_dp) then
      call random_number(u)
      x = random_gamma(shape + 1.0_dp) * max(u, tiny(1.0_dp))**(1.0_dp/shape)
      return
    end if

    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      z = random_normal()
      v = 1.0_dp + c*z
      if (v <= 0.0_dp) cycle
      v = v*v*v
      call random_number(u)
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(max(u, tiny(1.0_dp))) < 0.5_dp*z*z + d*(1.0_dp - v + log(v))) exit
    end do
    x = d*v
  end function random_gamma

  function random_chisq(df) result(x)
    real(dp), intent(in) :: df
    real(dp) :: x

    x = 2.0_dp*random_gamma(0.5_dp*df)
  end function random_chisq

  function random_student_t(df) result(x)
    real(dp), intent(in) :: df
    real(dp) :: x

    x = random_normal()/sqrt(random_chisq(df)/df)
  end function random_student_t

end module betategarch_rng
