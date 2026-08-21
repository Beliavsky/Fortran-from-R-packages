! SPDX-License-Identifier: GPL-2.0-only
module hypergeo_kinds
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    private
    integer, parameter, public :: dp = real64
    real(dp), parameter, public :: pi = acos(-1.0_dp)
    complex(dp), parameter, public :: ci = (0.0_dp, 1.0_dp)
end module hypergeo_kinds
