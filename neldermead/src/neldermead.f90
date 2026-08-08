! SPDX-License-Identifier: CECILL-2.0
! Derived from the R package neldermead 1.0-13 and its Scilab lineage.
! See LICENSE and UPSTREAM_PROVENANCE.md.

module neldermead
  use neldermead_kinds, only : dp
  use neldermead_types
  use neldermead_simplex
  use neldermead_core
  use neldermead_frontends
  implicit none
  public
end module neldermead
