! Modern Fortran translation of the computational code in the R package
! pmultinom 1.0.0 by Alexander Davis.
! Upstream license: GNU Affero General Public License v3 (AGPL-3).
module pmultinom_kinds
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)
end module pmultinom_kinds
