! SPDX-License-Identifier: BSD-3-Clause
program time_series_metrics
    use metrics, only : dp, mase, smape
    implicit none

    real(dp), parameter :: actual(6) = [1.1_dp, 1.9_dp, 3.0_dp, 4.4_dp, 5.0_dp, 5.6_dp]
    real(dp), parameter :: predicted(6) = [0.9_dp, 1.8_dp, 2.5_dp, 4.5_dp, 5.0_dp, 6.2_dp]

    write(*, '(a, f10.6)') 'MASE(1): ', mase(actual, predicted, 1)
    write(*, '(a, f10.6)') 'MASE(2): ', mase(actual, predicted, 2)
    write(*, '(a, f10.6)') 'SMAPE:   ', smape(actual, predicted)
end program time_series_metrics
