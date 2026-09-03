program test_fcs
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_fcs_result, mice_fcs_impute, mice_method_pmm, mice_method_norm, mice_method_logreg
    implicit none
    real(dp) :: data(10, 3), original(10, 3), nan
    integer :: methods(3), pred(3, 3), info, i, chain
    type(mice_fcs_result) :: result

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    do i = 1, 10
        data(i, 1) = real(i, dp)
        data(i, 2) = 2.0_dp + 0.5_dp * real(i, dp)
        if (i <= 5) then
            data(i, 3) = 0.0_dp
        else
            data(i, 3) = 1.0_dp
        end if
    end do
    data(3, 1) = nan
    data(7, 1) = nan
    data(2, 2) = nan
    data(8, 2) = nan
    data(4, 3) = nan
    data(9, 3) = nan
    original = data
    methods = [mice_method_pmm, mice_method_norm, mice_method_logreg]
    pred = 1
    do i = 1, 3
        pred(i, i) = 0
    end do
    call mice_fcs_impute(data, methods, pred, 3, 2, 2026_int64, result, info)
    if (info /= mice_ok) error stop "FCS status"
    if (any(ieee_is_nan(result%completed))) error stop "FCS left NaNs"
    do chain = 1, 2
        do i = 1, 10
            if (.not. ieee_is_nan(original(i, 1))) then
                if (abs(result%completed(i, 1, chain) - original(i, 1)) > 0.0_dp) error stop "FCS changed observed data"
            end if
        end do
    end do
    if (any(result%completed(:, 3, :) < 0.0_dp) .or. any(result%completed(:, 3, :) > 1.0_dp)) error stop "FCS binary support"
    print *, "test_fcs: PASS"
end program test_fcs
