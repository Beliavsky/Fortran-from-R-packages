! SPDX-License-Identifier: BSD-3-Clause
program demo_metrics
    use metrics
    implicit none

    real(dp), parameter :: actual_reg(4) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: predicted_reg(4) = [1.1_dp, 1.8_dp, 3.2_dp, 3.7_dp]
    integer, parameter :: actual_bin(4) = [1, 1, 0, 0]
    real(dp), parameter :: probability(4) = [0.9_dp, 0.7_dp, 0.3_dp, 0.1_dp]

    write(*, '(a)') 'Metrics modern Fortran demonstration'
    write(*, '(a, f10.6)') 'RMSE: ', rmse(actual_reg, predicted_reg)
    write(*, '(a, f10.6)') 'MAE:  ', mae(actual_reg, predicted_reg)
    write(*, '(a, f10.6)') 'AUC:  ', auc(actual_bin, probability)
    write(*, '(a, f10.6)') 'MAP@3:', apk(3, [1, 3, 5], [1, 2, 3, 4, 5])
end program demo_metrics
