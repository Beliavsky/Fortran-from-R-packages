! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
module manifoldoptim_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
end module manifoldoptim_kinds
