! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint_kinds
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter, public :: dp = real64
end module classint_kinds
