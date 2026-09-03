program test_fcs_parity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_fcs_result, mice_fcs_impute, mice_method_midastouch, mice_method_polr, &
                     mice_method_lda
    implicit none

    real(dp) :: data(30, 3), nan
    integer :: methods(3), predictor(3, 3), categories(3), info, i
    type(mice_fcs_result) :: result

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    do i = 1, 30
        data(i, 1) = 0.2_dp * real(i, dp) + sin(real(i, dp)) / 20.0_dp
        if (i <= 10) then
            data(i, 2) = 1.0_dp
            data(i, 3) = 1.0_dp
        else if (i <= 20) then
            data(i, 2) = 2.0_dp
            data(i, 3) = 2.0_dp
        else
            data(i, 2) = 3.0_dp
            data(i, 3) = 3.0_dp
        end if
    end do
    data(5, 1) = nan
    data(16, 1) = nan
    data(27, 1) = nan
    data(4, 2) = nan
    data(15, 2) = nan
    data(26, 2) = nan
    data(3, 3) = nan
    data(14, 3) = nan
    data(25, 3) = nan
    methods = [mice_method_midastouch, mice_method_polr, mice_method_lda]
    predictor = 1
    do i = 1, 3
        predictor(i, i) = 0
    end do
    categories = [0, 3, 3]
    call mice_fcs_impute(data, methods, predictor, 1, 1, 808_int64, result, info, category_count=categories)
    if (info /= mice_ok) error stop "parity FCS status"
    if (any(ieee_is_nan(result%completed))) error stop "parity FCS left NaN"
    if (any(result%completed(:, 2, 1) < 1.0_dp) .or. any(result%completed(:, 2, 1) > 3.0_dp)) then
        error stop "parity FCS polr support"
    end if
    if (any(result%completed(:, 3, 1) < 1.0_dp) .or. any(result%completed(:, 3, 1) > 3.0_dp)) then
        error stop "parity FCS lda support"
    end if
    print *, "test_fcs_parity: PASS"
end program test_fcs_parity
