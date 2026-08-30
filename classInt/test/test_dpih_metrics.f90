program test_dpih_metrics
    use classint
    implicit none

    type(class_intervals) :: fit
    type(classint_options) :: opts
    type(jenks_indices) :: indices
    real(dp) :: x(20)
    real(dp) :: area(20)
    real(dp) :: h
    integer :: i

    do i = 1, 20
        x(i) = real(i, dp) + 0.1_dp * sin(real(i, dp))
        area(i) = 1.0_dp + 0.05_dp * real(i, dp)
    end do
    h = dpih_bandwidth(x, "minim", 2, 401, [minval(x), maxval(x)], .true.)
    if (.not. (h > 0.0_dp .and. h < maxval(x) - minval(x))) error stop "dpih bandwidth"

    opts = classint_options()
    opts%dpih_level = 2
    call classint_fit(x, "dpih", fit, opts)
    if (size(fit%breaks) < 2) error stop "dpih class intervals"

    call classint_fit(x, 4, "fisher", fit)
    indices = jenks_tests(fit, area)
    if (.not. (indices%goodness_of_fit > 0.0_dp .and. indices%goodness_of_fit <= 1.0_dp)) &
        error stop "gvf range"
    if (.not. indices%has_overview) error stop "overview flag"
    if (.not. (classint_aic(fit) > -huge(1.0_dp))) error stop "finite AIC"

    print *, "test_dpih_metrics: ok"
end program test_dpih_metrics
