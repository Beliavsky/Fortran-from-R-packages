program test_core
    use classint
    implicit none

    type(class_intervals) :: fit
    type(classint_options) :: opts
    integer, allocatable :: cols(:)
    integer, allocatable :: counts(:)
    real(dp) :: x(6)
    real(dp) :: y(30)
    real(dp) :: ll
    real(dp), allocatable :: fish_stats(:, :)
    real(dp), allocatable :: fish_breaks(:)
    real(dp) :: fish_x(8)
    integer :: i

    x = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 50.0_dp]
    call classint_fit(x, 3, "fisher", fit)
    call assert_close_vec(fit%breaks, [0.0_dp, 0.5_dp, 26.0_dp, 50.0_dp], 1.0e-12_dp, "fisher breaks")
    if (size(fit%fisher_stats, 1) /= 3 .or. size(fit%fisher_stats, 2) /= 4) error stop "fisher stats shape"

    call classint_fit(x, 3, "jenks", fit)
    call assert_close_vec(fit%breaks, [0.0_dp, 0.0_dp, 2.0_dp, 50.0_dp], 1.0e-12_dp, "jenks breaks")
    if (trim(fit%interval_closure) /= "right") error stop "jenks must force right closure"
    call find_cols(fit, cols)
    if (any(cols /= [1, 1, 1, 2, 2, 3])) error stop "jenks class labels"

    call classify_intervals(x, 3, "jenks", cols)
    if (any(cols /= [1, 1, 1, 2, 2, 3])) error stop "fit-and-classify convenience"

    do i = 1, 10
        y(i) = 1.0_dp
        y(i + 10) = 2.0_dp
        y(i + 20) = 3.0_dp
    end do
    call classint_fit(y, 2, "jenks", fit)
    ll = classint_loglik(fit)
    if (abs(ll - (-14.52876_dp)) > 1.0e-5_dp) error stop "upstream logLik fixture"
    if (abs(classint_n_partitions(fit) - 2.0_dp) > 1.0e-12_dp) error stop "partition count fixture"

    call classint_fit(y, 3, "jenks", fit)
    if (abs(classint_loglik(fit)) > 1.0e-12_dp) error stop "singleton-within-class logLik fixture"

    fish_x = [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, 8.0_dp, 9.0_dp, 10.0_dp]
    call fisher_exact(fish_x, 3, fish_stats, fish_breaks)
    call assert_close_vec(fish_breaks, [1.0_dp, 3.0_dp, 6.5_dp, 10.0_dp], 1.0e-13_dp, &
                          "upstream fish break fixture")
    call assert_close_vec(fish_stats(1, :), [8.0_dp, 10.0_dp, 9.0_dp, 0.8164965809277289_dp], 5.0e-14_dp, &
                          "upstream fish high class")
    call assert_close_vec(fish_stats(2, :), [4.0_dp, 5.0_dp, 4.5_dp, 0.5_dp], 5.0e-14_dp, &
                          "upstream fish middle class")
    call assert_close_vec(fish_stats(3, :), [1.0_dp, 2.0_dp, 1.6666666666666667_dp, 0.4714045207910313_dp], &
                          5.0e-14_dp, "upstream fish low class")

    opts = classint_options()
    allocate(opts%fixed_breaks(4))
    opts%fixed_breaks = [-1.0_dp, 0.5_dp, 10.0_dp, 60.0_dp]
    call classint_fit(x, 3, "fixed", fit, opts)
    counts = class_counts(fit)
    if (any(counts /= [3, 2, 1])) error stop "fixed class counts"

    print *, "test_core: ok"
contains
    subroutine assert_close_vec(actual, expected, tol, label)
        real(dp), intent(in) :: actual(:) !! Computed vector checked elementwise against a deterministic regression fixture.
        real(dp), intent(in) :: expected(:) !! Expected regression vector of the same length as actual.
        real(dp), intent(in) :: tol !! Maximum allowed absolute elementwise difference.
        character(len=*), intent(in) :: label !! Short failure label emitted if the vector comparison fails.

        if (size(actual) /= size(expected)) error stop "assert_close_vec: size mismatch"
        if (maxval(abs(actual - expected)) > tol) then
            print *, trim(label), actual, expected
            error stop "assert_close_vec: mismatch"
        end if
    end subroutine assert_close_vec
end program test_core
