! SPDX-License-Identifier: GPL-3.0-only
module qpoases_kinds
    implicit none
    integer, parameter, public :: dp = kind(1.0d0)
    real(dp), parameter, public :: qpoases_infinity = huge(1.0_dp) / 100.0_dp
end module qpoases_kinds
