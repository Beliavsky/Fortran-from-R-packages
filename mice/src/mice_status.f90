! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Computational translation derived from mice 3.19.0.
module mice_status
    implicit none
    private

    integer, parameter, public :: mice_ok = 0
    integer, parameter, public :: mice_invalid_argument = -1
    integer, parameter, public :: mice_invalid_shape = -2
    integer, parameter, public :: mice_singular = 1
    integer, parameter, public :: mice_no_observed = 2
    integer, parameter, public :: mice_not_converged = 3

end module mice_status
