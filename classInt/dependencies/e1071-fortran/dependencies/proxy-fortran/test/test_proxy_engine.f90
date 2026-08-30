program test_proxy_engine
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use proxy, only: dp, proxy_dist_auto, proxy_dist_cross, proxy_dist_pairwise, &
                     proxy_simil_auto, proxy_measure_value, proxy_pack_dist, proxy_unpack_dist, &
                     proxy_subset_dist, proxy_row_sums_dist, proxy_row_means_dist, proxy_row_dist, proxy_col_dist
    implicit none

    real(dp) :: x(3, 2)
    real(dp) :: y(3, 2)
    real(dp), allocatable :: d(:, :)
    real(dp), allocatable :: s(:, :)
    real(dp), allocatable :: cross(:, :)
    real(dp), allocatable :: pair(:)
    real(dp), allocatable :: packed(:)
    real(dp), allocatable :: unpacked(:, :)
    real(dp), allocatable :: subset(:)
    real(dp), allocatable :: sums(:)
    real(dp), allocatable :: means(:)
    integer, allocatable :: rows(:)
    integer, allocatable :: cols(:)
    integer :: status
    real(dp) :: unknown_value

    x(1, :) = [0.0_dp, 0.0_dp]
    x(2, :) = [3.0_dp, 4.0_dp]
    x(3, :) = [3.0_dp, 0.0_dp]
    y = x
    y(:, 1) = y(:, 1) + 1.0_dp

    call proxy_dist_auto(x, 'Euclidean', d, status=status)
    call assert_true(status == 0, 'euclidean auto status')
    call assert_close(d(1, 2), 5.0_dp, 1.0e-12_dp, 'euclidean auto value')
    call assert_close(d(2, 3), 4.0_dp, 1.0e-12_dp, 'euclidean auto second value')

    call proxy_dist_cross(x, y, 'Manhattan', cross, status=status)
    call assert_true(status == 0, 'cross status')
    call assert_close(cross(1, 1), 1.0_dp, 1.0e-12_dp, 'cross manhattan')

    call proxy_dist_pairwise(x, y, 'Euclidean', pair, status=status)
    call assert_true(status == 0, 'pairwise status')
    call assert_true(all(abs(pair - 1.0_dp) < 1.0e-12_dp), 'pairwise euclidean')

    call proxy_simil_auto(x, 'Euclidean', s, status=status)
    call assert_close(s(1, 2), 1.0_dp / 6.0_dp, 1.0e-12_dp, 'distance to similarity')

    call proxy_dist_auto(merge(1.0_dp, 0.0_dp, x > 0.0_dp), 'Jaccard', s, status=status)
    call assert_close(s(2, 3), 0.5_dp, 1.0e-12_dp, 'jaccard converted to distance')

    call proxy_pack_dist(d, packed)
    call assert_true(size(packed) == 3, 'packed length')
    call assert_close(packed(1), 5.0_dp, 1.0e-12_dp, 'packed order 1')
    call assert_close(packed(2), 3.0_dp, 1.0e-12_dp, 'packed order 2')
    call assert_close(packed(3), 4.0_dp, 1.0e-12_dp, 'packed order 3')
    call proxy_unpack_dist(packed, 3, unpacked)
    call assert_true(maxval(abs(unpacked - d)) < 1.0e-12_dp, 'unpack parity')

    call proxy_subset_dist(packed, 3, [1, 3], subset)
    call assert_close(subset(1), 3.0_dp, 1.0e-12_dp, 'subset dist')
    call proxy_row_sums_dist(packed, 3, sums)
    call assert_true(maxval(abs(sums - [8.0_dp, 9.0_dp, 7.0_dp])) < 1.0e-12_dp, 'row sums')
    call proxy_row_means_dist(packed, 3, means, diag=.false.)
    call assert_true(maxval(abs(means - [4.0_dp, 4.5_dp, 3.5_dp])) < 1.0e-12_dp, 'row means')
    call proxy_row_dist(3, rows)
    call proxy_col_dist(3, cols)
    call assert_true(all(rows == [2, 3, 3]), 'row.dist indices')
    call assert_true(all(cols == [1, 1, 2]), 'col.dist indices')

    status = 0
    unknown_value = proxy_measure_value(x(1, :), x(2, :), 'not-a-method', .true., status=status)
    call assert_true(status == 1, 'unknown method status')
    call assert_true(ieee_is_nan(unknown_value), 'unknown method returns NaN')

contains

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! Calculated scalar result for this engine assertion.
        real(dp), intent(in) :: expected !! Deterministic expected value for the selected proxy operation.
        real(dp), intent(in) :: tolerance !! Allowed absolute difference between calculated and expected values.
        character(len=*), intent(in) :: label !! Descriptive assertion label printed on failure.

        if (abs(actual - expected) > tolerance) then
            write (*, '(a,2es24.14)') 'FAIL '//trim(label)//': ', actual, expected
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition !! Boolean condition that must hold for the test to pass.
        character(len=*), intent(in) :: label !! Descriptive assertion label printed when `condition` is false.

        if (.not. condition) then
            write (*, '(a)') 'FAIL '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

end program test_proxy_engine
