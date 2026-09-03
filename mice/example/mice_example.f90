program mice_example
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_fcs_result, mice_fcs_impute, mice_method_pmm, mice_method_norm
    implicit none
    real(dp) :: data(6, 2), nan
    integer :: methods(2), predictors(2, 2), info, i
    type(mice_fcs_result) :: result

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    data(:, 1) = [1.0_dp, 2.0_dp, nan, 4.0_dp, 5.0_dp, 6.0_dp]
    data(:, 2) = [2.0_dp, nan, 6.0_dp, 8.0_dp, 10.0_dp, 12.0_dp]
    methods = [mice_method_pmm, mice_method_norm]
    predictors = reshape([0, 1, 1, 0], [2, 2])
    call mice_fcs_impute(data, methods, predictors, 5, 2, 1234_int64, result, info)
    if (info /= mice_ok) error stop "mice example failed"
    print '(a)', "First completed data set:"
    do i = 1, size(data, 1)
        print '(2f12.5)', result%completed(i, :, 1)
    end do
end program mice_example
