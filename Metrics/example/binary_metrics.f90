! SPDX-License-Identifier: BSD-3-Clause
program binary_metrics
    use metrics, only : dp, auc, logloss, precision, recall, fbeta_score
    implicit none

    integer, parameter :: actual(6) = [1, 1, 1, 0, 0, 0]
    real(dp), parameter :: probability(6) = [0.9_dp, 0.8_dp, 0.4_dp, 0.5_dp, 0.3_dp, 0.2_dp]
    integer, parameter :: predicted(6) = [1, 1, 0, 1, 0, 0]

    write(*, '(a, f10.6)') 'AUC:       ', auc(actual, probability)
    write(*, '(a, f10.6)') 'Log loss:  ', logloss(actual, probability)
    write(*, '(a, f10.6)') 'Precision: ', precision(actual, predicted)
    write(*, '(a, f10.6)') 'Recall:    ', recall(actual, predicted)
    write(*, '(a, f10.6)') 'F1:        ', fbeta_score(actual, predicted)
end program binary_metrics
