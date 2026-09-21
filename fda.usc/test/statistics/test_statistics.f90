program test_statistics
    use fdausc
    implicit none

    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp) :: t(5)
    real(dp) :: x(4, 5)
    real(dp) :: y(4)
    real(dp) :: d(3, 3)
    real(dp) :: rpstat(2, 2)
    real(dp) :: residuals(2, 3)
    real(dp) :: inprod(3)
    real(dp) :: avec(4)
    integer :: order(3, 2)
    integer :: groups(4)
    real(dp) :: fanova

    inprod = [1.0_dp, 0.25_dp, 1.0_dp]
    call adot_matrix_vector(2, inprod, avec)
    if (any(avec < 0.0_dp)) error stop "Adot produced a negative packed weight"

    residuals(1, :) = [-1.0_dp, 0.5_dp, 0.5_dp]
    residuals(2, :) = [0.5_dp, 0.5_dp, -1.0_dp]
    order(:, 1) = [1, 2, 3]
    order(:, 2) = [3, 2, 1]
    call rp_projection_statistics(order, residuals, rpstat)
    if (any(rpstat < 0.0_dp)) error stop "random-projection statistics must be nonnegative"

    d = reshape([0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 0.0_dp], [3, 3])
    call assert_close(distance_correlation(d, d), 1.0_dp, tol, "distance correlation self check")

    t = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
    x(1, :) = t
    x(2, :) = t + 1.0_dp
    x(3, :) = 2.0_dp*t
    x(4, :) = 2.0_dp*t + 1.0_dp
    y = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    if (flm_f_statistic(t, x, y) <= 0.0_dp) error stop "FLM statistic should be positive"

    groups = [1, 1, 2, 2]
    fanova = fanova_mean_distance_statistic(t, x, groups)
    if (fanova <= 0.0_dp) error stop "functional ANOVA mean statistic should be positive"

    print *, "All fda.usc statistics tests passed"
contains
    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! Computed numerical value being checked.
        real(dp), intent(in) :: expected !! Reference value required by the deterministic test.
        real(dp), intent(in) :: tolerance !! Maximum permitted absolute error for the check.
        character(len=*), intent(in) :: label !! Short diagnostic label printed on failure.

        if (abs(actual - expected) > tolerance) then
            print *, trim(label), actual, expected
            error stop "numerical assertion failed"
        end if
    end subroutine assert_close
end program test_statistics
