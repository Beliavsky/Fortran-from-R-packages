! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

module rspectra_types
    use rspectra_kinds, only: dp
    implicit none
    private

    type, public :: eigs_opts
        integer :: ncv = 0
        integer :: maxitr = 1000
        real(dp) :: tol = 1.0e-10_dp
        logical :: retvec = .true.
        real(dp), allocatable :: initvec(:)
    end type eigs_opts

    type, public :: eigs_sym_result
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: vectors(:,:)
        integer :: nconv = 0
        integer :: niter = 0
        integer :: nops = 0
        integer :: info = 0
    end type eigs_sym_result

    type, public :: eigs_result
        complex(dp), allocatable :: values(:)
        complex(dp), allocatable :: vectors(:,:)
        integer :: nconv = 0
        integer :: niter = 0
        integer :: nops = 0
        integer :: info = 0
    end type eigs_result

    type, public :: svds_opts
        integer :: ncv = 0
        integer :: maxitr = 1000
        real(dp) :: tol = 1.0e-10_dp
    end type svds_opts

    type, public :: svds_result
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: u(:,:)
        real(dp), allocatable :: v(:,:)
        integer :: nconv = 0
        integer :: niter = 0
        integer :: nops = 0
        integer :: info = 0
    end type svds_result

end module rspectra_types
