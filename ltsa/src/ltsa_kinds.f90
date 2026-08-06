! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_kinds
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)
    real(dp), parameter, public :: pi = acos(-1.0_dp)
    real(dp), parameter, public :: two_pi = 2.0_dp*pi
end module ltsa_kinds
