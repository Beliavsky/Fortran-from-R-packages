! SPDX-License-Identifier: BSD-3-Clause
program test_time_series
    use metrics, only : dp, mase
    use test_support
    implicit none

    real(dp), parameter :: actual(6) = [1.1_dp, 1.9_dp, 3.0_dp, 4.4_dp, 5.0_dp, 5.6_dp]
    real(dp), parameter :: predicted(6) = [0.9_dp, 1.8_dp, 2.5_dp, 4.5_dp, 5.0_dp, 6.2_dp]

    call check_close('mase lag 1', mase(actual, predicted, 1), 1.5_dp / 5.4_dp)
    call check_close('mase lag 2', mase(actual, predicted, 2), 1.5_dp / 11.4_dp)

    call finish_tests('test_time_series')
end program test_time_series
