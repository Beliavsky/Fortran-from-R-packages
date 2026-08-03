! SPDX-License-Identifier: BSD-3-Clause
program test_regression
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use metrics, only : dp, bias, percent_bias, se, sse, mse, rmse, ae, mae, mdae, ape, mape, smape, &
                        sle, msle, rmsle, rae, rse, rrse, explained_variation
    use test_support
    implicit none

    real(dp) :: nan_value, pos_inf, neg_inf

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)
    neg_inf = ieee_value(0.0_dp, ieee_negative_inf)

    call check_close('bias scalar', bias([1.0_dp], [1.0_dp]), 0.0_dp)
    call check_close('bias vector', bias([-1.0_dp, -100.0_dp, 17.5_dp], [0.0_dp, 0.0_dp, 0.0_dp]), -27.833333333333333_dp)
    call check_close('percent bias', percent_bias([1.0_dp, 2.0_dp, 3.0_dp], [1.0_dp, 3.0_dp, 2.0_dp]), &
                     (-0.5_dp + 1.0_dp / 3.0_dp) / 3.0_dp)
    call check_close('percent bias negative infinity', percent_bias([1.0_dp, 2.0_dp, 0.0_dp], [1.0_dp, 2.0_dp, 1.0_dp]), neg_inf)
    call check_close('percent bias nan', percent_bias([0.0_dp], [0.0_dp]), nan_value)
    call check_array_close('squared error', se([9.0_dp, 10.0_dp, 11.0_dp], [11.0_dp, 10.0_dp, 9.0_dp]), [4.0_dp, 0.0_dp, 4.0_dp])
    call check_close('sse', sse([1.0_dp, 3.0_dp, 2.0_dp], [2.0_dp, 3.0_dp, 4.0_dp]), 5.0_dp)
    call check_close('mse', mse([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2.0_dp, 3.0_dp, 4.0_dp, 4.0_dp]), 0.75_dp)
    call check_close('rmse', rmse([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp]), 0.5_dp)
    call check_array_close('absolute error', ae([9.0_dp, 10.0_dp, 11.0_dp], [11.0_dp, 10.0_dp, 9.0_dp]), [2.0_dp, 0.0_dp, 2.0_dp])
    call check_close('mae', mae([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp]), 0.25_dp)
    call check_close('mdae', mdae([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 4.0_dp, 50.0_dp]), 0.5_dp)
    call check_array_close('ape', ape([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]), &
                           [pos_inf, 1.0_dp, 0.5_dp, 1.0_dp / 3.0_dp])
    call check_array_close('ape nan', ape([0.0_dp, 1.0_dp, 2.0_dp], [0.0_dp, 0.0_dp, 0.0_dp]), [nan_value, 1.0_dp, 1.0_dp])
    call check_close('mape', mape([1.0_dp, 2.0_dp, 3.0_dp], [2.0_dp, 3.0_dp, 4.0_dp]), (1.0_dp + 0.5_dp + 1.0_dp / 3.0_dp) / 3.0_dp)
    call check_close('smape nan', smape([0.0_dp], [0.0_dp]), nan_value)
    call check_close('smape opposite', smape([1.0_dp], [-1.0_dp]), 2.0_dp)
    call check_array_close('sle', sle([0.0_dp, 1.0_dp, 3.4_dp], [1.0_dp, 0.0_dp, 3.4_dp]), &
                           [log(2.0_dp)**2, log(2.0_dp)**2, 0.0_dp])
    call check_close('msle', msle([1.0_dp, 2.0_dp, exp(1.0_dp) - 1.0_dp], &
                                  [1.0_dp, 2.0_dp, exp(2.0_dp) - 1.0_dp]), 1.0_dp / 3.0_dp)
    call check_close('rmsle', rmsle([exp(5.0_dp) - 1.0_dp], [exp(1.0_dp) - 1.0_dp]), 4.0_dp)
    call check_close('rae', rae([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp]), 0.25_dp)
    call check_close('rse', rse([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp]), 0.2_dp)
    call check_close('rrse', rrse([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp]), sqrt(0.2_dp))
    call check_close('explained variation', explained_variation([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
                                                                [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp]), 0.8_dp)

    call finish_tests('test_regression')
end program test_regression
