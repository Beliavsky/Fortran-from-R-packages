! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of fNonlinear/tseriesChaos computational routines and is distributed
! under the GNU General Public License version 2 or later.
module chaos_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
end module chaos_kinds
