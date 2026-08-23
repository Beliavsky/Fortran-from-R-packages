! SPDX-License-Identifier: GPL-3.0-only
module bigstatsr_kinds
    use iso_fortran_env, only: real64, int8, int32, int64
    implicit none
    private
    integer, parameter, public :: dp = kind(1.0d0)
    public :: int8, int32, int64
end module bigstatsr_kinds
