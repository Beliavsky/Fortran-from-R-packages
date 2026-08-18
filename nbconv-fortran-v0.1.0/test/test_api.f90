! SPDX-License-Identifier: GPL-3.0-or-later
program test_api
    use nbconv, only : dp, dnbconv, dnbconv_mu, dnbconv_p, pnbconv_p, qnbconv_p
    implicit none
    real(dp), parameter :: ps(2) = [0.05_dp, 0.44_dp]
    real(dp), parameter :: phis(2) = [5.0_dp, 8.0_dp]
    real(dp), parameter :: mus(2) = [95.0_dp, 10.1818181818181818_dp]
    integer, parameter :: counts(3) = [0, 50, 100]
    integer, parameter :: qs(4) = [25, 50, 100, 200]
    real(dp), parameter :: cdf_ref(4) = [ &
        3.53712535669463771e-03_dp, 7.45977690078236666e-02_dp, &
        5.15605412583506673e-01_dp, 9.68208253760752702e-01_dp]
    real(dp), parameter :: probs(5) = [0.05_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.95_dp]
    integer, parameter :: qref(5) = [46, 73, 99, 130, 186]
    real(dp), allocatable :: a(:), b(:), c(:), g(:)
    integer, allocatable :: q(:)

    a = dnbconv_mu(counts, mus, phis, method="exact")
    b = dnbconv_p(counts, ps, phis, method="exact")
    call assert_close_vec(a, b, 2.0e-13_dp, 1.0e-15_dp)
    g = dnbconv(counts, mus, phis, parameterization="mu", method="exact")
    call assert_close_vec(a, g, 2.0e-13_dp, 1.0e-15_dp)

    c = pnbconv_p(qs, ps, phis, method="exact")
    call assert_close_vec(c, cdf_ref, 3.0e-12_dp, 2.0e-14_dp)

    q = qnbconv_p(probs, 500, ps, phis, method="exact")
    if (any(q /= qref)) error stop 1
    print *, "test_api: PASS"

contains

    subroutine assert_close_vec(actual, expected, rtol, atol)
        real(dp), intent(in) :: actual(:), expected(:), rtol, atol
        integer :: i
        if (size(actual) /= size(expected)) error stop 1
        do i = 1, size(actual)
            if (abs(actual(i) - expected(i)) > atol + rtol * abs(expected(i))) error stop 1
        end do
    end subroutine assert_close_vec

end program test_api
