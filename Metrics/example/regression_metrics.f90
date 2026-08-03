! SPDX-License-Identifier: BSD-3-Clause
program regression_metrics
    use metrics, only : dp, mae, rmse, mape, explained_variation
    implicit none

    real(dp), parameter :: actual(6) = [1.1_dp, 1.9_dp, 3.0_dp, 4.4_dp, 5.0_dp, 5.6_dp]
    real(dp), parameter :: predicted(6) = [0.9_dp, 1.8_dp, 2.5_dp, 4.5_dp, 5.0_dp, 6.2_dp]

    write(*, '(a, f10.6)') 'MAE:  ', mae(actual, predicted)
    write(*, '(a, f10.6)') 'RMSE: ', rmse(actual, predicted)
    write(*, '(a, f10.6)') 'MAPE: ', mape(actual, predicted)
    write(*, '(a, f10.6)') 'R2:   ', explained_variation(actual, predicted)
end program regression_metrics
