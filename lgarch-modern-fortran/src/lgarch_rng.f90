! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
module lgarch_rng
  use lgarch_kinds, only : dp, pi
  implicit none
  private
  public :: seed_rng, random_normal, random_normal_vector
contains
  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: values(:)
    call random_seed(size=n)
    allocate(values(n))
    do i=1,n
      values(i) = modulo(seed + 104729*i, huge(1)-1)
      if (values(i) <= 0) values(i) = i
    end do
    call random_seed(put=values)
  end subroutine seed_rng

  real(dp) function random_normal() result(z)
    real(dp) :: u1, u2
    call random_number(u1); call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function random_normal

  subroutine random_normal_vector(z)
    real(dp), intent(out) :: z(:)
    integer :: i
    do i=1,size(z)
      z(i) = random_normal()
    end do
  end subroutine random_normal_vector
end module lgarch_rng
