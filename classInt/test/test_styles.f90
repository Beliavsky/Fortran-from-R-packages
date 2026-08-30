program test_styles
    use classint
    implicit none

    type(class_intervals) :: fit
    type(class_intervals) :: fit2
    type(classint_options) :: opts
    real(dp) :: x(12)
    real(dp) :: heavy(10)
    integer, allocatable :: cols(:)

    x = [1.0_dp, 1.2_dp, 1.3_dp, 2.0_dp, 2.1_dp, 2.2_dp, &
         8.0_dp, 8.1_dp, 8.2_dp, 12.0_dp, 12.1_dp, 12.2_dp]

    call classint_fit(x, 3, "equal", fit)
    if (size(fit%breaks) /= 4) error stop "equal breaks"
    call classint_fit(x, 3, "pretty", fit)
    if (size(fit%breaks) < 3) error stop "pretty breaks"
    call classint_fit(x, 3, "quantile", fit)
    if (abs(fit%breaks(2) - 2.0_dp) > 0.2_dp) error stop "quantile break"

    opts = classint_options()
    opts%seed = 17
    opts%kmeans_nstart = 5
    call classint_fit(x, 3, "kmeans", fit, opts)
    call find_cols(fit, cols)
    if (count(cols == 1) /= 6 .or. count(cols == 2) /= 3 .or. count(cols == 3) /= 3) &
        error stop "kmeans groups"

    opts = classint_options()
    opts%hclust_method = "complete"
    call classint_fit(x, 3, "hclust", fit, opts)
    if (.not. fit%has_hclust) error stop "hclust model missing"
    call get_hclust_class_intervals(fit, 4, fit2)
    if (size(fit2%breaks) /= 5) error stop "hclust recut"

    opts = classint_options()
    opts%seed = 9
    opts%bclust_iter_base = 3
    opts%bclust_base_centers = 3
    opts%bclust_maxcluster = 3
    opts%bclust_resample = .false.
    opts%bclust_hclust = "average"
    call classint_fit(x, 3, "bclust", fit, opts)
    if (.not. fit%has_bclust) error stop "bclust model missing"
    if (size(fit%breaks) /= 4) error stop "bclust breaks"

    call classint_fit(x, 3, "sd", fit)
    if (size(fit%breaks) < 3) error stop "sd breaks"

    heavy = [1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 8.0_dp, 13.0_dp, 34.0_dp]
    call classint_fit(heavy, "headtails", fit)
    if (fit%breaks(1) > minval(heavy) .or. fit%breaks(size(fit%breaks)) < maxval(heavy)) &
        error stop "headtails endpoints"

    call classint_fit(x, 3, "maximum", fit)
    if (size(fit%breaks) < 4) error stop "maximum breaks"

    opts = classint_options()
    opts%box_iqr_mult = 1.5_dp
    call classint_fit(x, "box", fit, opts)
    if (size(fit%breaks) /= 7) error stop "box breaks"

    print *, "test_styles: ok"
end program test_styles
