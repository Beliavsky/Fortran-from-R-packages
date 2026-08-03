! SPDX-License-Identifier: BSD-3-Clause
program test_binary_classification
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use metrics, only : dp, auc, ll, logloss, precision, recall, fbeta_score
    use test_support
    implicit none

    real(dp) :: inf
    inf = ieee_value(0.0_dp, ieee_positive_inf)

    call check_close('auc 1', auc([1, 0, 1, 1], [0.32_dp, 0.52_dp, 0.26_dp, 0.86_dp]), 1.0_dp / 3.0_dp)
    call check_close('auc perfect', auc([1, 0, 1, 0, 1], [0.9_dp, 0.1_dp, 0.8_dp, 0.1_dp, 0.7_dp]), 1.0_dp)
    call check_close('auc ties', auc([1, 1, 1, 1, 0, 0, 0, 0, 0, 0], &
        [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]), 0.5_dp)
    call check_close('ll correct one', ll(1, 1.0_dp), 0.0_dp)
    call check_close('ll impossible one', ll(1, 0.0_dp), inf)
    call check_close('ll impossible zero', ll(0, 1.0_dp), inf)
    call check_close('ll half', ll(1, 0.5_dp), -log(0.5_dp))
    call check_close('logloss exact', logloss([1, 1, 0, 0], [1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]), 0.0_dp)
    call check_close('logloss reference', logloss([1, 1, 1, 0, 0, 0], &
        [0.5_dp, 0.1_dp, 0.01_dp, 0.9_dp, 0.75_dp, 0.001_dp]), 1.881797068998267_dp, 1.0e-14_dp)
    call check_close('precision', precision([1, 1, 0, 0], [1, 1, 1, 1]), 0.5_dp)
    call check_close('recall', recall([1, 1, 1, 1], [1, 0, 0, 1]), 0.5_dp)
    call check_close('fbeta', fbeta_score([0, 0, 1, 1], [1, 1, 1, 0]), 0.4_dp)
    call check_close('fbeta beta zero', fbeta_score([1, 1, 0, 0], [1, 1, 1, 1], beta=0.0_dp), 0.5_dp)

    call finish_tests('test_binary_classification')
end program test_binary_classification
