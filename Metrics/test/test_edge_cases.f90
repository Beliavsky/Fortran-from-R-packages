! SPDX-License-Identifier: BSD-3-Clause
program test_edge_cases
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use metrics
    use test_support
    implicit none

    integer :: stat
    real(dp) :: value

    value = mse([1.0_dp], [1.0_dp, 2.0_dp], stat)
    call check_true('mismatched regression is nan', ieee_is_nan(value))
    call check_true('mismatched regression status', stat == metrics_invalid_size)

    value = auc([1, 1], [0.1_dp, 0.2_dp])
    call check_true('auc one class is nan', ieee_is_nan(value))

    value = precision([1, 0], [0, 0])
    call check_true('precision no predicted positives is nan', ieee_is_nan(value))

    value = apk(0, [1], [1])
    call check_true('apk zero k is nan', ieee_is_nan(value))

    value = mase([1.0_dp, 2.0_dp], [1.0_dp, 2.0_dp], 2, stat)
    call check_true('invalid mase is nan', ieee_is_nan(value))
    call check_true('invalid mase status', stat == metrics_invalid_argument)

    call finish_tests('test_edge_cases')
end program test_edge_cases
