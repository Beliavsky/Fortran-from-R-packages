! SPDX-License-Identifier: GPL-3.0-or-later
program test_degenerate
    use nbconv, only : dp, nb_sum_exact, dnbconv_mu
    implicit none
    real(dp), parameter :: ps(2) = [1.0_dp, 1.0_dp]
    real(dp), parameter :: phis(2) = [2.0_dp, 7.0_dp]
    real(dp), parameter :: mus(2) = [0.0_dp, 0.0_dp]
    integer, parameter :: counts(4) = [0, 1, 2, 10]
    real(dp), allocatable :: a(:), b(:)

    a = nb_sum_exact(ps, phis, counts)
    b = dnbconv_mu(counts, mus, phis, method="moments")
    if (abs(a(1) - 1.0_dp) > 1.0e-15_dp .or. any(a(2:) > 0.0_dp)) error stop 1
    if (abs(b(1) - 1.0_dp) > 1.0e-15_dp .or. any(b(2:) > 0.0_dp)) error stop 1
    print *, "test_degenerate: PASS"
end program test_degenerate
