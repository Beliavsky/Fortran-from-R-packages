! Modern Fortran translation of computational code from TSP 1.2.7.
! Original Copyright (C) Michael Hahsler and Kurt Hornik.
! SPDX-License-Identifier: GPL-3.0-only
! See LICENSE, COPYING, and UPSTREAM.md for provenance and licensing.

module tsp_kinds
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)

end module tsp_kinds
