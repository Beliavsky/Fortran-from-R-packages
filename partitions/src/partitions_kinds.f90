! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

module partitions_kinds
    use, intrinsic :: iso_fortran_env, only : int64, real64
    implicit none
    private

    integer, parameter, public :: i8 = int64
    integer, parameter, public :: dp = real64

end module partitions_kinds
