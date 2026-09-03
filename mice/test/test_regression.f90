program test_regression
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_rng_state, rng_seed, estimice, impute_norm, impute_norm_predict
    implicit none
    real(dp) :: x(8, 1), design(8, 2), y(8)
    real(dp), allocatable :: coef(:), residuals(:), vcov(:, :), imp(:), pred(:)
    logical :: observed(8), where(8), ridge_used
    type(mice_rng_state) :: rng
    integer :: df, i, info

    do i = 1, 8
        x(i, 1) = real(i - 1, dp)
        design(i, :) = [1.0_dp, x(i, 1)]
        y(i) = 2.0_dp + 3.0_dp * x(i, 1)
    end do
    call estimice(design, y, 1.0e-5_dp, coef, residuals, vcov, df, ridge_used, info)
    if (info /= mice_ok) error stop "estimice status"
    if (maxval(abs(coef - [2.0_dp, 3.0_dp])) > 1.0e-11_dp) error stop "estimice coefficients"
    if (maxval(abs(residuals)) > 1.0e-11_dp) error stop "estimice residuals"
    observed = [.true., .true., .true., .true., .true., .true., .false., .false.]
    where = .not. observed
    call rng_seed(rng, 9876_int64)
    call impute_norm(y, observed, x, where, rng, imp, info)
    if (info /= mice_ok .or. size(imp) /= 2) error stop "normal imputation"
    call impute_norm_predict(y, observed, x, where, pred, info)
    if (info /= mice_ok) error stop "predict status"
    if (maxval(abs(pred - [20.0_dp, 23.0_dp])) > 1.0e-9_dp) error stop "normal prediction"
    print *, "test_regression: PASS"
end program test_regression
