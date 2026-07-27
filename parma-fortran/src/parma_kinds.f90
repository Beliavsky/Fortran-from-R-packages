! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Derived from parma 1.7, Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
module parma_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: parma_small = 1.0e-20_dp
   real(dp), parameter, public :: parma_big = 1.0e20_dp
   real(dp), parameter, public :: pi = acos(-1.0_dp)
end module parma_kinds
