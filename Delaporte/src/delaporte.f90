! Copyright (c) 2013, Avraham Adler
! All rights reserved.
! SPDX-License-Identifier: BSD-2-Clause
!
! Public facade for the standalone Delaporte library.

module delaporte
    use delaporte_kinds, only : dp
    use delaporte_distribution, only : ddelap, pdelap, qdelap, ddelap_vec, &
        pdelap_vec, qdelap_vec
    use delaporte_rng, only : rdelap, rdelap_vec, seed_delaporte
    use delaporte_moments, only : momdelap
    use delaporte_approx, only : qdelap_approx
    implicit none
    private

    public :: dp
    public :: ddelap, pdelap, qdelap
    public :: ddelap_vec, pdelap_vec, qdelap_vec
    public :: rdelap, rdelap_vec, seed_delaporte
    public :: momdelap, qdelap_approx

end module delaporte
