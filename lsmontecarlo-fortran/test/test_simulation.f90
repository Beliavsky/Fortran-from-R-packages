! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
program test_simulation
    use lsmc_kinds, only : dp
    use lsmc_random, only : seed_random_number
    use lsmc_simulation, only : fast_gbm, first_value_row
    use lsmc_simulation, only : simulate_antithetic_correlated_gbm_paths
    use lsmc_simulation, only : simulate_antithetic_gbm_paths
    implicit none

    real(dp), allocatable :: paths(:, :)
    real(dp) :: pair1(1000, 12)
    real(dp) :: pair2(1000, 12)
    real(dp) :: single_pair(1000, 12)
    real(dp) :: expected_product
    real(dp) :: x(3, 4)
    real(dp) :: y(3, 4)
    integer :: j
    integer :: positive_counts(3)

    paths = fast_gbm(spot=100.0_dp, sigma=0.0_dp, n=4, m=10, rate=0.03_dp, dividend=0.01_dp, &
        maturity=1.0_dp, seed=123)
    call assert_true(all(shape(paths) == [4, 10]), "fast_gbm shape")
    call assert_true(all(paths > 0.0_dp), "fast_gbm positivity")
    call assert_close(paths(1, 10), 100.0_dp * exp(0.02_dp), 2.0e-12_dp, "deterministic terminal GBM")
    call assert_true(maxval(abs(paths - spread(paths(1, :), 1, 4))) < 2.0e-12_dp, "deterministic paths")

    call seed_random_number(456)
    call simulate_antithetic_gbm_paths(100.0_dp, 0.2_dp, 500, 12, 0.04_dp, 0.01_dp, 1.0_dp, single_pair)
    do j = 1, 12
        expected_product = 10000.0_dp * exp(2.0_dp * (0.04_dp - 0.01_dp - 0.5_dp * 0.2_dp**2) * &
            real(j, dp) / 12.0_dp)
        call assert_true(maxval(abs(single_pair(1:500, j) * single_pair(501:1000, j) - expected_product)) < 1.0e-8_dp, &
            "antithetic GBM identity")
    end do

    call seed_random_number(789)
    call simulate_antithetic_correlated_gbm_paths(100.0_dp, 0.2_dp, 0.04_dp, 0.01_dp, &
        100.0_dp, 0.2_dp, 0.04_dp, 0.01_dp, 1.0_dp, 500, 12, 1.0_dp, pair1, pair2)
    call assert_true(maxval(abs(pair1 - pair2)) < 1.0e-12_dp, "perfectly correlated paths")

    x = reshape([0.0_dp, 0.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, &
        0.0_dp, 7.0_dp, 8.0_dp], [3, 4])
    y = first_value_row(x)
    positive_counts = count(y > 0.0_dp, dim=2)
    call assert_true(all(positive_counts == 1), "one retained value per row")
    call assert_close(sum(y), 10.0_dp, 0.0_dp, "first value row values")

    print '(a)', 'test_simulation: PASS'

contains

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            print '(a)', trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tolerance
        character(len=*), intent(in) :: label

        call assert_true(abs(actual - expected) <= tolerance, label)
    end subroutine assert_close

end program test_simulation
