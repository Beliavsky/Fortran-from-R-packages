! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

module timsac_raw
  use timsac_kinds, only: dp
  implicit none
  private

  public :: auspecf
  public :: autcorf
  public :: mulcorf
  public :: wnoisef

  interface
    subroutine auspecf(n, lagh1, cxx1, p1, p2, q)
      import :: dp
      implicit none
      integer :: n, lagh1
      real(dp) :: cxx1(lagh1), p1(lagh1), p2(lagh1), q(lagh1)
    end subroutine auspecf

    subroutine autcorf(x, n, cxx, cn, lagh1, xmean)
      import :: dp
      implicit none
      integer :: n, lagh1
      real(dp) :: x(n), cxx(lagh1), cn(lagh1), xmean
    end subroutine autcorf

    subroutine mulcorf(x1, n, k, lagh1, sm, c, cn)
      import :: dp
      implicit none
      integer :: n, k, lagh1
      real(dp) :: x1(n,k), sm(k), c(lagh1,k,k), cn(lagh1,k,k)
    end subroutine mulcorf

    subroutine wnoisef(nra, ir, sd1, x2)
      import :: dp
      implicit none
      integer :: nra, ir
      real(dp) :: sd1(ir,ir), x2(ir,nra)
    end subroutine wnoisef
  end interface

end module timsac_raw
