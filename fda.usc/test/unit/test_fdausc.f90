program test_fdausc
    use fdausc
    implicit none

    real(dp), parameter :: tol = 2.0e-8_dp
    real(dp) :: x(5)
    real(dp) :: y(5)
    real(dp) :: curves(3, 5)
    real(dp) :: mean_curve(5)
    real(dp) :: variance_curve(5)
    real(dp) :: distances(3, 3)
    real(dp) :: s(3, 3)
    real(dp) :: pvalues(3)
    real(dp) :: measures(3)
    real(dp) :: components(1, 5)
    real(dp) :: scores(3, 1)
    real(dp) :: singular_values(1)
    real(dp) :: pca_mean(5)
    real(dp) :: phi(3, 2)
    real(dp) :: penalty(2, 2)
    real(dp) :: sbasis(3, 3)
    integer :: info
    integer :: i

    x = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
    y = x*x
    call assert_close(integrate_curve(x, y, 2), 1.0_dp/3.0_dp, 2.0e-4_dp, "Simpson integration")
    call assert_close(ker_epa(0.0_dp), 0.75_dp, tol, "Epanechnikov kernel")
    call assert_close(iker_unif(0.0_dp), 0.5_dp, tol, "integrated uniform kernel")

    curves(1, :) = x
    curves(2, :) = 2.0_dp*x
    curves(3, :) = 3.0_dp*x
    call functional_mean(curves, mean_curve)
    call functional_variance(curves, variance_curve)
    call assert_close(mean_curve(5), 2.0_dp, tol, "functional mean")
    call assert_close(variance_curve(5), 2.0_dp/3.0_dp, tol, "functional variance")

    call lp_distance_matrix(x, curves, curves, distances, p=2.0_dp, method=2)
    call assert_close(distances(1, 1), 0.0_dp, tol, "Lp distance diagonal")
    call assert_close(distances(1, 3), sqrt(4.0_dp/3.0_dp), 2.0e-3_dp, "Lp distance")

    call smoothing_nw([0.0_dp, 1.0_dp, 2.0_dp], 1.0_dp, s, 1)
    do i = 1, 3
        call assert_close(sum(s(i, :)), 1.0_dp, 5.0e-12_dp, "NW row sum")
    end do

    pvalues = [0.001_dp, 0.02_dp, 0.8_dp]
    if (.not. fdr_reject(pvalues, 0.95_dp, 1)) error stop "FDR rejection test failed"

    call prediction_measures([1.0_dp, 2.0_dp, 4.0_dp], [1.0_dp, 3.0_dp, 2.0_dp], measures)
    call assert_close(measures(1), sqrt(5.0_dp/3.0_dp), tol, "prediction RMSE")
    call assert_close(measures(2), 1.0_dp, tol, "prediction MAE")
    call assert_close(measures(3), 5.0_dp/3.0_dp, tol, "prediction MSE")

    call functional_pca(x, curves, 1, components, scores, singular_values, pca_mean, info)
    if (info /= 0) error stop "functional PCA decomposition failed"
    if (singular_values(1) <= 0.0_dp) error stop "functional PCA singular value should be positive"
    call assert_close(sum(scores(:, 1)), 0.0_dp, 5.0e-10_dp, "centered PCA scores")

    phi(:, 1) = 1.0_dp
    phi(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp]
    penalty = 0.0_dp
    call basis_smoothing_matrix(phi, [1.0_dp, 1.0_dp, 1.0_dp], penalty, 0.0_dp, sbasis, info)
    if (info /= 0) error stop "basis smoother solve failed"
    call assert_close(matrix_trace(sbasis), 2.0_dp, 5.0e-10_dp, "basis smoother degrees of freedom")

    print *, "All fda.usc unit tests passed"
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
end program test_fdausc
