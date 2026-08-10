! SPDX-License-Identifier: GPL-2.0-or-later
module bvls_types
    use bvls_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: bvls_success = 0
    integer, parameter, public :: bvls_max_iterations = 1
    integer, parameter, public :: bvls_invalid_dimensions = -1
    integer, parameter, public :: bvls_inconsistent_bounds = -2
    integer, parameter, public :: bvls_no_free_variables = -3
    integer, parameter, public :: bvls_invalid_state = -4
    integer, parameter, public :: bvls_rank_failure = -5

    type, public :: bvls_result
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: fitted(:)
        real(dp), allocatable :: residuals(:)
        real(dp), allocatable :: gradient(:)
        integer, allocatable :: istate(:)
        real(dp) :: deviance = 0.0_dp
        real(dp) :: residual_norm = 0.0_dp
        integer :: iterations = 0
        integer :: status = bvls_success
    end type bvls_result

end module bvls_types
