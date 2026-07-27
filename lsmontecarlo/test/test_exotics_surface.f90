! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
program test_exotics_surface
    use lsmontecarlo, only : amer_put_lsm, amer_put_lsm_price_surface
    use lsmontecarlo, only : asian_amer_put_lsm, dp, option_result, price_surface
    use lsmontecarlo, only : quanto_amer_put_lsm, quanto_amer_put_lsm_av
    use lsmontecarlo, only : surface_maximum, surface_mean, surface_minimum
    implicit none

    type(option_result) :: asian
    type(option_result) :: plain
    type(option_result) :: quanto
    type(option_result) :: quanto_av
    type(price_surface) :: surface
    real(dp) :: strikes(2)
    real(dp) :: volatilities(2)

    plain = amer_put_lsm(100.0_dp, 0.25_dp, 25000, 40, 105.0_dp, 0.04_dp, 0.01_dp, 1.0_dp, 2024)
    asian = asian_amer_put_lsm(100.0_dp, 0.25_dp, 25000, 40, 105.0_dp, 0.04_dp, 0.01_dp, 1.0_dp, 2024)
    quanto = quanto_amer_put_lsm(100.0_dp, 0.25_dp, 25000, 40, 105.0_dp, 0.04_dp, 0.01_dp, 1.0_dp, &
        2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2024)
    quanto_av = quanto_amer_put_lsm_av(100.0_dp, 0.25_dp, 12500, 40, 105.0_dp, 0.04_dp, 0.01_dp, 1.0_dp, &
        2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2024)

    call assert_true(asian%price > 0.0_dp .and. asian%price < 105.0_dp, "Asian price bounds")
    call assert_true(abs(quanto%price - 2.0_dp * plain%price) < 0.45_dp, "constant multiplier quanto")
    call assert_true(abs(quanto_av%price - 2.0_dp * plain%price) < 0.45_dp, "antithetic constant multiplier quanto")
    call assert_true(quanto%standard_error > 0.0_dp, "quanto standard error")

    volatilities = [0.15_dp, 0.30_dp]
    strikes = [95.0_dp, 110.0_dp]
    surface = amer_put_lsm_price_surface(volatilities, strikes, spot=100.0_dp, n=4000, m=25, &
        rate=0.04_dp, dividend=0.01_dp, maturity=1.0_dp, seed=404)
    call assert_true(all(shape(surface%values) == [2, 2]), "surface shape")
    call assert_true(all(surface%values >= 0.0_dp), "surface nonnegative")
    call assert_true(surface_minimum(surface) <= surface_mean(surface), "surface minimum")
    call assert_true(surface_mean(surface) <= surface_maximum(surface), "surface maximum")
    call assert_true(maxval(abs(surface%volatilities - volatilities)) <= epsilon(1.0_dp), "surface volatility labels")
    call assert_true(maxval(abs(surface%strikes - strikes)) <= epsilon(1.0_dp), "surface strike labels")

    print '(a)', 'test_exotics_surface: PASS'

contains

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            print '(a)', trim(label)
            print '(a,4f12.6)', 'prices: ', plain%price, asian%price, quanto%price, quanto_av%price
            error stop 1
        end if
    end subroutine assert_true

end program test_exotics_surface
