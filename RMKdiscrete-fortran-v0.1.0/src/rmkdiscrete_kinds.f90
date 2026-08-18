! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
module rmkdiscrete_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
end module rmkdiscrete_kinds
