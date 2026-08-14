! FNN-fortran: modern Fortran translation of computational code from FNN 1.1.4.1.
! Modified/translated 2026 by the FNN-fortran contributors.
! SPDX-License-Identifier: GPL-2.0-or-later
! See UPSTREAM.md and upstream/FNN-1.1.4.1 for original authorship and notices.
module fnn_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
end module fnn_kinds
