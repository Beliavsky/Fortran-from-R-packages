module test_proxy_callbacks
    use proxy, only: dp
    implicit none
    private
    public :: scaled_l1

contains

    function scaled_l1(a, b, p) result(value)
        real(dp), intent(in) :: a(:) !! First numeric observation supplied by the custom-registry test.
        real(dp), intent(in) :: b(:) !! Second numeric observation supplied by the custom-registry test.
        real(dp), intent(in), optional :: p(:) !! Optional scale vector; this fixture uses `p(1)` when present.
        real(dp) :: value
        real(dp) :: scale

        scale = 1.0_dp
        if (present(p)) then
            if (size(p) > 0) scale = p(1)
        end if
        value = scale * sum(abs(a - b))
    end function scaled_l1

end module test_proxy_callbacks

program test_proxy_special
    use proxy, only: dp, binary_similarity_from_counts, nominal_similarity, &
                     gower_auto_similarity, proxy_gower_logical, proxy_gower_factor, &
                     proxy_gower_metric, proxy_gower_ordinal, proxy_register_numeric, &
                     proxy_registry_clear, proxy_dist_auto, proxy_numeric_callback
    use test_proxy_callbacks, only: scaled_l1
    implicit none

    real(dp) :: gx(3, 4)
    integer :: types(4)
    real(dp), allocatable :: g(:, :)
    real(dp) :: nx(6)
    real(dp) :: ny(6)
    real(dp) :: x(2, 2)
    real(dp), allocatable :: d(:, :)
    integer :: status
    procedure(proxy_numeric_callback), pointer :: custom_pointer

    call assert_close(binary_similarity_from_counts('Jaccard', 2, 1, 1, 2, 6), &
                      0.5_dp, 1.0e-12_dp, 'binary jaccard')
    call assert_close(binary_similarity_from_counts('Jaccard', 0, 0, 0, 6, 6), &
                      1.0_dp, 1.0e-12_dp, 'binary jaccard all false')
    call assert_close(binary_similarity_from_counts('Dice', 2, 1, 1, 2, 6), &
                      2.0_dp / 3.0_dp, 1.0e-12_dp, 'binary dice')
    call assert_close(binary_similarity_from_counts('Simple Matching', 2, 1, 1, 2, 6), &
                      2.0_dp / 3.0_dp, 1.0e-12_dp, 'simple matching')

    nx = [1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 3.0_dp]
    ny = [1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp]
    call assert_close(nominal_similarity(nx, ny, 'Phi-squared'), &
                      37.0_dp / 36.0_dp, 1.0e-12_dp, 'phi squared')

    gx(1, :) = [1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]
    gx(2, :) = [0.0_dp, 1.0_dp, 10.0_dp, 2.0_dp]
    gx(3, :) = [1.0_dp, 2.0_dp, 5.0_dp, 3.0_dp]
    types = [proxy_gower_logical, proxy_gower_factor, proxy_gower_metric, proxy_gower_ordinal]
    call gower_auto_similarity(gx, types, g)
    call assert_close(g(1, 1), 1.0_dp, 1.0e-12_dp, 'gower diagonal')
    call assert_true(g(1, 2) < 1.0_dp .and. g(1, 2) >= 0.0_dp, 'gower range')
    call assert_true(maxval(abs(g - transpose(g))) < 1.0e-12_dp, 'gower symmetry')

    call proxy_registry_clear()
    custom_pointer => scaled_l1
    call proxy_register_numeric('scaled-l1', custom_pointer, .true., status=status)
    call assert_true(status == 0, 'registry registration')
    x(1, :) = [0.0_dp, 0.0_dp]
    x(2, :) = [2.0_dp, 3.0_dp]
    call proxy_dist_auto(x, 'scaled-l1', d, p=[2.0_dp], status=status)
    call assert_true(status == 0, 'registry dispatch')
    call assert_close(d(1, 2), 10.0_dp, 1.0e-12_dp, 'custom registry value')

contains

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! Calculated scalar result for this specialized-measure assertion.
        real(dp), intent(in) :: expected !! Reference value expected from the upstream measure definition.
        real(dp), intent(in) :: tolerance !! Maximum allowed absolute numerical error.
        character(len=*), intent(in) :: label !! Descriptive assertion label emitted on failure.

        if (abs(actual - expected) > tolerance) then
            write (*, '(a,2es24.14)') 'FAIL '//trim(label)//': ', actual, expected
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition !! Boolean test condition that must evaluate to true.
        character(len=*), intent(in) :: label !! Human-readable label emitted when the assertion fails.

        if (.not. condition) then
            write (*, '(a)') 'FAIL '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

end program test_proxy_special
