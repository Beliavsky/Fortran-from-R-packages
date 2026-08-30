program classint_example
    use classint
    implicit none

    type(class_intervals) :: fit
    real(dp) :: x(10)
    integer, allocatable :: classes(:)
    integer :: i

    x = [2.0_dp, 3.0_dp, 3.5_dp, 4.0_dp, 8.0_dp, 9.0_dp, 10.0_dp, 17.0_dp, 18.0_dp, 30.0_dp]
    call classint_fit(x, 4, "fisher", fit)
    call classify_intervals(fit, classes)

    print '(a,*(f8.3,1x))', "breaks: ", fit%breaks
    do i = 1, size(x)
        print '(f8.3,2x,i0)', x(i), classes(i)
    end do
end program classint_example
