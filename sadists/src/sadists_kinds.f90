! SPDX-License-Identifier: LGPL-3.0-or-later
module sadists_kinds
    implicit none
    private
    integer, parameter, public :: dp = kind(1.0d0)
    real(dp), parameter, public :: pi = acos(-1.0_dp)
    real(dp), parameter, public :: sqrt2 = sqrt(2.0_dp)
    real(dp), parameter, public :: log2 = log(2.0_dp)
end module sadists_kinds
