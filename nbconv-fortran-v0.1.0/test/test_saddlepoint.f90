! SPDX-License-Identifier: GPL-3.0-or-later
program test_saddlepoint
    use nbconv, only : dp, nb_sum_saddlepoint
    implicit none
    real(dp), parameter :: mus(2) = [95.0_dp, 10.1818181818181818_dp]
    real(dp), parameter :: phis(2) = [5.0_dp, 8.0_dp]
    integer :: counts(0:500), i
    integer, parameter :: idx(8) = [0, 1, 10, 25, 50, 100, 200, 500]
    real(dp), parameter :: ref(8) = [ &
        4.33001129603604630e-10_dp, 4.33508075511694559e-09_dp, &
        1.15446473760885254e-05_dp, 6.49304157201850999e-04_dp, &
        5.37649667458264553e-03_dp, 9.39590431067686450e-03_dp, &
        1.06515752062205237e-03_dp, 9.62060266158549380e-09_dp]
    real(dp), allocatable :: pmf(:)

    counts = [(i, i=0,500)]
    pmf = nb_sum_saddlepoint(mus, phis, counts, normalize=.true.)
    call assert_close(sum(pmf), 1.0_dp, 1.0e-13_dp, 1.0e-13_dp)
    do i = 1, size(idx)
        call assert_close(pmf(idx(i) + 1), ref(i), 3.0e-8_dp, 1.0e-12_dp)
    end do
    print *, "test_saddlepoint: PASS"

contains

    subroutine assert_close(actual, expected, rtol, atol)
        real(dp), intent(in) :: actual, expected, rtol, atol
        if (abs(actual - expected) > atol + rtol * abs(expected)) error stop 1
    end subroutine assert_close

end program test_saddlepoint
