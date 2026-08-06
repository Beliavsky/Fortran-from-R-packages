! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
module matlab_types
    use matlab_kinds, only : dp
    implicit none
    private

    type, public :: meshgrid2d_result
        real(dp), allocatable :: x(:, :)
        real(dp), allocatable :: y(:, :)
    end type meshgrid2d_result

    type, public :: meshgrid3d_result
        real(dp), allocatable :: x(:, :, :)
        real(dp), allocatable :: y(:, :, :)
        real(dp), allocatable :: z(:, :, :)
    end type meshgrid3d_result

    type, public :: fileparts_result
        character(len=:), allocatable :: pathstr
        character(len=:), allocatable :: name
        character(len=:), allocatable :: ext
        character(len=:), allocatable :: versn
    end type fileparts_result
end module matlab_types
