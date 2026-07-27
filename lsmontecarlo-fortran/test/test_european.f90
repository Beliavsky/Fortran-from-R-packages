! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
program test_european
    use lsmontecarlo, only : dp, eu_call_bs, eu_put_bs
    implicit none

    real(dp) :: call_price
    real(dp) :: parity_error
    real(dp) :: put_price

    call_price = eu_call_bs(100.0_dp, 0.20_dp, 100.0_dp, 0.05_dp, 0.0_dp, 1.0_dp)
    put_price = eu_put_bs(100.0_dp, 0.20_dp, 100.0_dp, 0.05_dp, 0.0_dp, 1.0_dp)
    call assert_close(call_price, 10.450583572185565_dp, 2.0e-12_dp, "Black-Scholes call")
    call assert_close(put_price, 5.573526022256971_dp, 2.0e-12_dp, "Black-Scholes put")

    parity_error = call_price - put_price - (100.0_dp - 100.0_dp * exp(-0.05_dp))
    call assert_close(parity_error, 0.0_dp, 2.0e-12_dp, "put-call parity")
    call assert_close(eu_call_bs(100.0_dp, 0.0_dp, 90.0_dp, 0.05_dp, 0.0_dp, 1.0_dp), &
        exp(-0.05_dp) * max(100.0_dp * exp(0.05_dp) - 90.0_dp, 0.0_dp), 1.0e-12_dp, "zero-volatility call")
    call assert_close(eu_put_bs(90.0_dp, 0.20_dp, 100.0_dp, 0.05_dp, 0.0_dp, 0.0_dp), &
        10.0_dp, 0.0_dp, "maturity payoff")

    print '(a)', 'test_european: PASS'

contains

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tolerance
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tolerance) then
            print '(a,2(1x,es24.16))', trim(label), actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_european
