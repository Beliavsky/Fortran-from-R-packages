! SPDX-License-Identifier: GPL-2.0-only
module hypergeo_types
    use hypergeo_kinds, only : dp
    implicit none
    private

    type, public :: hypergeo_info
        logical :: converged = .false.
        integer :: iterations = 0
        integer :: method_code = 0
        real(dp) :: residual = huge(1.0_dp)
        character(len=32) :: method = ''
    end type hypergeo_info
end module hypergeo_types
