! SPDX-License-Identifier: GPL-2.0-only
!
! Modern Fortran translation derived from CRAN signal 1.8-1.
module signal_kinds
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    private

    integer, parameter, public :: dp = real64
    real(dp), parameter, public :: signal_pi = acos(-1.0_dp)

end module signal_kinds
