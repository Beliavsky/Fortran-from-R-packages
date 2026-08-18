! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from R package good 1.0.2.

module good_kinds
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)

end module good_kinds
