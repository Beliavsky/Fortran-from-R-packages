! SPDX-License-Identifier: GPL-2.0-or-later
module tensora_kinds
    use, intrinsic :: iso_fortran_env, only : real64, int64
    implicit none
    private

    integer, parameter, public :: dp = real64
    integer, parameter, public :: i8 = int64
    integer, parameter, public :: axis_name_len = 64

end module tensora_kinds
