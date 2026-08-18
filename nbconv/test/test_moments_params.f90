! SPDX-License-Identifier: GPL-3.0-or-later
program test_moments_params
    use nbconv, only : dp, nb_sum_moments, nbconv_summary, nbconv_params_mu
    implicit none
    real(dp), parameter :: mus(2) = [95.0_dp, 10.1818181818181818_dp]
    real(dp), parameter :: phis(2) = [5.0_dp, 8.0_dp]
    integer, parameter :: counts(5) = [0, 25, 50, 100, 200]
    real(dp), parameter :: ref(5) = [ &
        2.0876224992356732e-08_dp, 8.436147137430339e-04_dp, &
        5.312136173229187e-03_dp, 9.356982396904835e-03_dp, &
        1.063250518135893e-03_dp]
    real(dp), allocatable :: pmf(:)
    type(nbconv_summary) :: s

    pmf = nb_sum_moments(mus, phis, counts)
    call assert_close_vec(pmf, ref, 2.0e-12_dp, 1.0e-14_dp)
    s = nbconv_params_mu(mus, phis)
    call assert_close(s%mean, 105.18181818181819_dp, 2.0e-14_dp, 2.0e-14_dp)
    call assert_close(s%variance, 1923.1404958677685_dp, 2.0e-14_dp, 2.0e-12_dp)
    call assert_close(s%skewness, 0.8795940313271249_dp, 2.0e-13_dp, 2.0e-14_dp)
    call assert_close(s%excess_kurtosis, 1.1719239291964865_dp, 2.0e-13_dp, 2.0e-14_dp)
    call assert_close(s%k_mean, 69.64285714285712_dp, 2.0e-13_dp, 2.0e-13_dp)
    print *, "test_moments_params: PASS"

contains

    subroutine assert_close(actual, expected, rtol, atol)
        real(dp), intent(in) :: actual, expected, rtol, atol
        if (abs(actual - expected) > atol + rtol * abs(expected)) error stop 1
    end subroutine assert_close

    subroutine assert_close_vec(actual, expected, rtol, atol)
        real(dp), intent(in) :: actual(:), expected(:), rtol, atol
        integer :: i
        if (size(actual) /= size(expected)) error stop 1
        do i = 1, size(actual)
            call assert_close(actual(i), expected(i), rtol, atol)
        end do
    end subroutine assert_close_vec

end program test_moments_params
