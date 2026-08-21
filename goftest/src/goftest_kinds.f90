! SPDX-License-Identifier: GPL-2.0-or-later
module goftest_kinds
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)
    real(dp), parameter, public :: pi = 3.141592653589793238462643383279502884197_dp

end module goftest_kinds
