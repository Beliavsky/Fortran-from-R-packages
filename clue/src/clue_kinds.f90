! SPDX-License-Identifier: GPL-2.0-only
module clue_kinds
    implicit none
    private
    integer, parameter, public :: dp = kind(1.0d0)
    real(dp), parameter, public :: clue_eps = epsilon(1.0_dp)
end module clue_kinds
