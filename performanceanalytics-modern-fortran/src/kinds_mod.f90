! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module kinds_mod
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
end module kinds_mod
