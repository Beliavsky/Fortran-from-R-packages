program test_desirability
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use desirability, only : dp, d_max, d_min, d_target, d_box, d_arb, &
        d_categorical, d_overall, d_max_type, d_min_type, d_target_type, &
        d_box_type, d_arb_type, d_categorical_type, d_overall_type, &
        hold, predict, predict_all, numeric_input, categorical_input, &
        missing_input, desirability_input, quiet_nan
    implicit none

    type(d_max_type) :: dmax_obj, dmax_tol
    type(d_min_type) :: dmin_obj
    type(d_target_type) :: dtarget_obj
    type(d_box_type) :: dbox_obj
    type(d_arb_type) :: darb_obj
    type(d_categorical_type) :: dcat_obj
    type(d_overall_type) :: overall, mixed_overall
    real(dp), allocatable :: y(:), oa(:, :), ov(:)
    type(desirability_input) :: mixed(2, 2)
    real(dp) :: xref(5)
    real(dp), parameter :: eps = 2.0e-12_dp

    dmax_obj = d_max(1.0_dp, 3.0_dp, scale=2.0_dp)
    y = predict(dmax_obj, [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
    call assert_vector(y, [0.0_dp, 0.0_dp, 0.25_dp, 1.0_dp, 1.0_dp], eps, "d_max")
    call assert_close(dmax_obj%missing, 328350.0_dp / 980100.0_dp, eps, &
        "d_max missing")

    dmin_obj = d_min(1.0_dp, 3.0_dp, scale=2.0_dp)
    y = predict(dmin_obj, [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
    call assert_vector(y, [1.0_dp, 1.0_dp, 0.25_dp, 0.0_dp, 0.0_dp], eps, "d_min")
    call assert_close(dmin_obj%missing, dmax_obj%missing, eps, "d_min missing")

    dtarget_obj = d_target(1.0_dp, 2.0_dp, 4.0_dp, &
        low_scale=2.0_dp, high_scale=1.0_dp)
    y = predict(dtarget_obj, [0.0_dp, 1.0_dp, 1.5_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp])
    call assert_vector(y, [0.0_dp, 0.0_dp, 0.25_dp, 1.0_dp, 0.5_dp, 0.0_dp, 0.0_dp], &
        eps, "d_target")

    dbox_obj = d_box(-1.0_dp, 1.0_dp)
    y = predict(dbox_obj, [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp])
    call assert_vector(y, [0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], eps, "d_box")
    call assert_close(dbox_obj%missing, 1.0_dp, eps, "d_box missing")

    darb_obj = d_arb([2.0_dp, 0.0_dp, 1.0_dp], [1.0_dp, 0.0_dp, 0.5_dp])
    y = predict(darb_obj, [-1.0_dp, 0.5_dp, 1.5_dp, 3.0_dp])
    call assert_vector(y, [0.0_dp, 0.25_dp, 0.75_dp, 1.0_dp], eps, "d_arb")
    call assert_close(darb_obj%missing, 0.5_dp, eps, "d_arb missing")

    dcat_obj = d_categorical([character(len=6) :: "value1", "value2", "value3"], &
        [0.1_dp, 0.9_dp, 0.2_dp])
    y = predict(dcat_obj, [character(len=6) :: "value2", "value1"])
    call assert_vector(y, [0.9_dp, 0.1_dp], eps, "d_categorical")
    call assert_close(dcat_obj%missing, 0.4_dp, eps, "d_categorical missing")

    dmax_tol = d_max(1.0_dp, 3.0_dp, tol=0.05_dp)
    call assert_close(predict(dmax_tol, 0.0_dp), 0.05_dp, eps, "tolerance")
    call assert_close(predict(dmax_tol, quiet_nan()), 0.5_dp, eps, "missing default")
    if (.not. ieee_is_nan(predict(dmax_tol, quiet_nan(), missing=quiet_nan()))) then
        error stop "missing override should preserve NaN"
    end if

    overall = d_overall([hold(d_max(80.0_dp, 97.0_dp)), &
        hold(d_target(55.0_dp, 57.5_dp, 60.0_dp))])
    xref = 0.0_dp
    ov = predict(overall, reshape([81.09_dp, 59.85_dp], [1, 2]))
    call assert_close(ov(1), sqrt((1.09_dp / 17.0_dp) * 0.06_dp), eps, &
        "overall numeric")
    oa = predict_all(overall, reshape([81.09_dp, 59.85_dp], [1, 2]))
    call assert_close(oa(1, 1), 1.09_dp / 17.0_dp, eps, "overall d1")
    call assert_close(oa(1, 2), 0.06_dp, eps, "overall d2")
    call assert_close(oa(1, 3), ov(1), eps, "overall combined")

    mixed_overall = d_overall()
    call mixed_overall%add(d_max(0.0_dp, 10.0_dp))
    call mixed_overall%add(dcat_obj)
    mixed(1, 1) = numeric_input(5.0_dp)
    mixed(1, 2) = categorical_input("value2")
    mixed(2, 1) = missing_input()
    mixed(2, 2) = categorical_input("")
    ov = predict(mixed_overall, mixed)
    call assert_close(ov(1), sqrt(0.5_dp * 0.9_dp), eps, "overall mixed")
    call assert_close(ov(2), sqrt(0.5_dp * 0.4_dp), eps, "overall mixed missing")

    print '(a)', "All desirability tests passed."

contains

    subroutine assert_close(actual, expected, tol, label)
        real(dp), intent(in) :: actual, expected, tol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tol) then
            print '(a,2(1x,es24.16))', trim(label), actual, expected
            error stop "assert_close failed"
        end if
    end subroutine assert_close

    subroutine assert_vector(actual, expected, tol, label)
        real(dp), intent(in) :: actual(:), expected(:), tol
        character(len=*), intent(in) :: label
        integer :: i

        if (size(actual) /= size(expected)) error stop "assert_vector size mismatch"
        do i = 1, size(actual)
            call assert_close(actual(i), expected(i), tol, label)
        end do
    end subroutine assert_vector

end program test_desirability
