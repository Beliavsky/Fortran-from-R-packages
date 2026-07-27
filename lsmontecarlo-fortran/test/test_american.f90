! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
program test_american
    use lsmontecarlo, only : amer_put_lsm, amer_put_lsm_av, amer_put_lsm_cv, dp
    use lsmontecarlo, only : eu_put_bs, option_result, price
    implicit none

    type(option_result) :: antithetic
    type(option_result) :: controlled
    type(option_result) :: plain
    real(dp) :: benchmark
    real(dp) :: european

    benchmark = american_put_binomial(100.0_dp, 100.0_dp, 0.05_dp, 0.0_dp, 0.20_dp, 1.0_dp, 1200)
    european = eu_put_bs(100.0_dp, 0.20_dp, 100.0_dp, 0.05_dp, 0.0_dp, 1.0_dp)
    plain = amer_put_lsm(100.0_dp, 0.20_dp, 30000, 50, 100.0_dp, 0.05_dp, 0.0_dp, 1.0_dp, 1001)
    antithetic = amer_put_lsm_av(100.0_dp, 0.20_dp, 15000, 50, 100.0_dp, 0.05_dp, 0.0_dp, 1.0_dp, 1001)
    controlled = amer_put_lsm_cv(100.0_dp, 0.20_dp, 30000, 50, 100.0_dp, 0.05_dp, 0.0_dp, 1.0_dp, 1001)

    call assert_true(abs(plain%price - benchmark) < 0.30_dp, "plain LSMC versus binomial")
    call assert_true(abs(antithetic%price - benchmark) < 0.30_dp, "antithetic LSMC versus binomial")
    call assert_true(abs(controlled%price - benchmark) < 0.30_dp, "control-variate LSMC versus binomial")
    call assert_true(plain%price >= european - 0.05_dp, "American value above European value")
    call assert_true(plain%price <= 100.0_dp, "American upper bound")
    call assert_true(plain%standard_error > 0.0_dp, "plain standard error")
    call assert_true(antithetic%effective_paths == 30000, "antithetic effective path count")
    call assert_true(antithetic%n_paths == 15000, "antithetic original path count")
    call assert_true(abs(price(controlled) - controlled%price) <= epsilon(1.0_dp), "price accessor")

    print '(a)', 'test_american: PASS'

contains

    function american_put_binomial(spot, strike, rate, dividend, sigma, maturity, n_steps) result(value)
        real(dp), intent(in) :: spot
        real(dp), intent(in) :: strike
        real(dp), intent(in) :: rate
        real(dp), intent(in) :: dividend
        real(dp), intent(in) :: sigma
        real(dp), intent(in) :: maturity
        integer, intent(in) :: n_steps
        real(dp), allocatable :: option(:)
        real(dp) :: discount
        real(dp) :: down
        real(dp) :: dt
        real(dp) :: probability
        real(dp) :: stock
        real(dp) :: up
        real(dp) :: value
        integer :: i
        integer :: j

        dt = maturity / real(n_steps, dp)
        up = exp(sigma * sqrt(dt))
        down = 1.0_dp / up
        discount = exp(-rate * dt)
        probability = (exp((rate - dividend) * dt) - down) / (up - down)
        allocate(option(0:n_steps))
        do j = 0, n_steps
            stock = spot * up**j * down**(n_steps - j)
            option(j) = max(strike - stock, 0.0_dp)
        end do
        do i = n_steps - 1, 0, -1
            do j = 0, i
                stock = spot * up**j * down**(i - j)
                option(j) = max(strike - stock, discount * &
                    (probability * option(j + 1) + (1.0_dp - probability) * option(j)))
            end do
        end do
        value = option(0)
    end function american_put_binomial

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            print '(a)', trim(label)
            print '(a,4f12.6)', 'values: ', plain%price, antithetic%price, controlled%price, benchmark
            error stop 1
        end if
    end subroutine assert_true

end program test_american
