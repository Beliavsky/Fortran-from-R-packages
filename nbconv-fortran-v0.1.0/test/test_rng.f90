! SPDX-License-Identifier: GPL-3.0-or-later
program test_rng
    use nbconv, only : dp, nbconv_seed, rnbconv_mu
    implicit none
    real(dp), parameter :: mus(2) = [95.0_dp, 10.1818181818181818_dp]
    real(dp), parameter :: phis(2) = [5.0_dp, 8.0_dp]
    integer, allocatable :: x(:)
    real(dp) :: m, v, expected_m, expected_v

    call nbconv_seed(20260816_8)
    x = rnbconv_mu(50000, mus, phis)
    if (any(x < 0)) error stop 1
    m = sum(real(x, dp)) / real(size(x), dp)
    v = sum((real(x, dp) - m)**2) / real(size(x) - 1, dp)
    expected_m = 105.18181818181819_dp
    expected_v = 1923.1404958677685_dp
    if (abs(m - expected_m) > 0.8_dp) error stop 1
    if (abs(v - expected_v) > 120.0_dp) error stop 1
    print *, "test_rng: PASS"
end program test_rng
