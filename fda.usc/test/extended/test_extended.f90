program test_extended
    use fdausc
    implicit none

    real(dp), parameter :: tol = 1.0e-9_dp
    real(dp) :: t(5)
    real(dp) :: curves(3, 5)
    real(dp) :: basis(2, 5)
    real(dp) :: coefs(3, 2)
    real(dp) :: mean_curve(5)
    real(dp) :: reconstructed(3, 5)
    real(dp) :: h(3)
    real(dp) :: dmat(3, 3)
    real(dp) :: fdist(3, 3)
    real(dp) :: simdepth(1)
    real(dp) :: point(1, 2)
    real(dp) :: triangle(3, 2)
    real(dp) :: ldepth(2)
    real(dp) :: ded(2, 3)
    real(dp) :: drr(3, 3)
    real(dp) :: probs(2, 2)
    real(dp) :: errs(2)
    real(dp) :: noise(3, 5, 2)
    real(dp) :: bootmeans(2, 5)
    real(dp) :: bootdist(2)
    real(dp) :: center(5)
    real(dp) :: dband
    real(dp) :: ydist(4, 4)
    real(dp) :: xcov(4, 3)
    real(dp) :: dc(3)
    real(dp) :: stat
    real(dp) :: df1
    real(dp) :: df2
    integer :: pred(2)
    integer :: truth(4)
    integer :: cand(4, 2)
    integer :: idx(3, 2)
    integer :: best
    integer :: info
    integer :: i
    logical :: sel(3)

    t = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
    curves(1, :) = 1.0_dp + 2.0_dp*t
    curves(2, :) = -1.0_dp + 0.5_dp*t
    curves(3, :) = 3.0_dp - t
    basis(1, :) = 1.0_dp
    basis(2, :) = t
    call basis_project_grid(curves, basis, coefs, mean_curve, center=.false., info=info)
    if (info /= 0) error stop "basis projection failed"
    call basis_reconstruct_grid(coefs, basis, 0.0_dp*t, reconstructed)
    call assert_close(maxval(abs(reconstructed - curves)), 0.0_dp, 2.0e-12_dp, "basis reconstruction")

    call lp_distance_matrix(t, curves, curves, dmat, p=2.0_dp, method=2)
    call default_bandwidth_grid(dmat, 0.25_dp, 0.75_dp, h)
    if (any(h <= 0.0_dp)) error stop "default bandwidths must be positive"

    call horizontal_shift_distance_matrix(curves, curves, t, fdist)
    call assert_close(maxval(abs([(fdist(i, i), i=1, 3)])), 0.0_dp, tol, "horizontal-shift diagonal")

    call fourier_semimetric(t, curves, curves, 3, 0, 1.0_dp, fdist, info)
    if (info /= 0) error stop "Fourier semimetric fit failed"
    call assert_close(maxval(abs([(fdist(i, i), i=1, 3)])), 0.0_dp, 5.0e-9_dp, "Fourier semimetric diagonal")

    point(1, :) = [0.25_dp, 0.25_dp]
    triangle(1, :) = [0.0_dp, 0.0_dp]
    triangle(2, :) = [1.0_dp, 0.0_dp]
    triangle(3, :) = [0.0_dp, 1.0_dp]
    call bivariate_simplicial_depth(point, triangle, simdepth)
    call assert_close(simdepth(1), 1.0_dp, tol, "simplicial depth inside triangle")

    ded(1, :) = [0.0_dp, 1.0_dp, 2.0_dp]
    ded(2, :) = [2.0_dp, 1.0_dp, 0.0_dp]
    drr = reshape([0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 0.0_dp], [3, 3])
    call likelihood_depth(ded, 1.0_dp, ldepth, scale=.true., distances_ref_ref=drr)
    if (any(ldepth <= 0.0_dp)) error stop "likelihood depth should be positive"

    call depth_classify(reshape([0.8_dp, 0.2_dp, 0.1_dp, 0.9_dp], [2, 2]), pred, probs)
    if (any(pred /= [1, 2])) error stop "maximum-depth classifier failed"

    truth = [1, 1, 2, 2]
    cand(:, 1) = [1, 2, 2, 2]
    cand(:, 2) = truth
    call select_classification_cv(truth, cand, best, errs)
    if (best /= 2) error stop "classification CV selection failed"
    call assert_close(errs(2), 0.0_dp, tol, "classification CV error")

    idx(:, 1) = [1, 2, 3]
    idx(:, 2) = [1, 1, 1]
    noise = 0.0_dp
    call bootstrap_mean_replicates(t, curves, idx, noise, 0.5_dp, center, bootmeans, bootdist, dband)
    call assert_close(maxval(abs(bootmeans(1, :) - center)), 0.0_dp, tol, "bootstrap identity mean")
    if (dband < 0.0_dp) error stop "bootstrap radius must be nonnegative"

    xcov(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    xcov(:, 2) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    xcov(:, 3) = [3.0_dp, 2.0_dp, 1.0_dp, 0.0_dp]
    do i = 1, 4
        ydist(i, :) = abs(xcov(i, 1) - xcov(:, 1))
    end do
    call local_distance_correlation_select(xcov, ydist, 0.1_dp, dc, sel)
    call assert_close(dc(1), 1.0_dp, 1.0e-10_dp, "LMDC distance correlation")
    if (.not. any(sel)) error stop "LMDC local maximum selection failed"

    call hetero_anova_onefactor([0.0_dp, 0.1_dp, 2.0_dp, 2.2_dp], [1, 1, 2, 2], stat, df1, df2)
    if (stat <= 0.0_dp .or. df1 <= 0.0_dp .or. df2 <= 0.0_dp) error stop "heteroscedastic ANOVA core failed"

    print *, "All fda.usc extended tests passed"
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
end program test_extended
