! rmutil computational translation
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: sqrt2 = sqrt(2.0_dp)
   real(dp), parameter, public :: log2pi = log(2.0_dp*pi)
end module rmutil_kinds
