! Copyright (c) 2013, Avraham Adler
! All rights reserved.
! SPDX-License-Identifier: BSD-2-Clause
!
! Standalone modern Fortran adaptation of the Delaporte R package.

module delaporte_kinds
    use, intrinsic :: iso_fortran_env, only : real64, int64
    implicit none
    private

    integer, parameter, public :: dp = real64
    integer, parameter, public :: i64 = int64

end module delaporte_kinds
