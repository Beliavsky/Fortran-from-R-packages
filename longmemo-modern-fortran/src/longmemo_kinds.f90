! Part of the modern Fortran translation of longmemo 1.1-4.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original longmemo authors retain copyright; see ORIGINAL_PACKAGE.txt.
! SPDX-License-Identifier: GPL-2.0-or-later

module longmemo_kinds
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)
    real(dp), parameter, public :: pi = acos(-1.0_dp)
    real(dp), parameter, public :: two_pi = 2.0_dp*pi

end module longmemo_kinds
