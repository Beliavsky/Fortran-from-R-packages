! genalg-fortran -- translation of the computational core of R package genalg.
! Original package: Copyright (C) Egon Willighagen, Michel Ballings, contributors.
! License: GPL-2.0-only. See LICENSE and original/ for provenance.
module genalg_kinds
    implicit none
    private
    integer, parameter, public :: dp = kind(1.0d0)
end module genalg_kinds
