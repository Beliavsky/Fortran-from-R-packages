! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
program test_array
    use matlab, only : dp, eye, fliplr, flipud, rot90, repmat, reshape2d, &
                       meshgrid, meshgrid2d_result, padarray, shape_of, size_dim, &
                       numel, ndims, isempty, find
    use test_support
    implicit none

    real(dp), allocatable :: a(:, :), b(:, :), expected(:, :), v(:)
    integer, allocatable :: idx(:), dims(:)
    type(meshgrid2d_result) :: grid

    a = eye(4, 2)
    call assert_true(all(shape(a) == [4, 2]), 'eye shape')
    call assert_close(a(1, 1) + a(2, 2), 2.0_dp, 0.0_dp, 'eye diagonal')

    deallocate(a)
    allocate(a(2, 3))
    a = reshape([1.0_dp, 4.0_dp, 2.0_dp, 5.0_dp, 3.0_dp, 6.0_dp], [2, 3])
    b = fliplr(a)
    expected = reshape([3.0_dp, 6.0_dp, 2.0_dp, 5.0_dp, 1.0_dp, 4.0_dp], [2, 3])
    call assert_all_close(b, expected, 0.0_dp, 'fliplr')
    b = flipud(a)
    expected = reshape([4.0_dp, 1.0_dp, 5.0_dp, 2.0_dp, 6.0_dp, 3.0_dp], [2, 3])
    call assert_all_close(b, expected, 0.0_dp, 'flipud')
    b = rot90(a)
    expected = reshape([3.0_dp, 2.0_dp, 1.0_dp, 6.0_dp, 5.0_dp, 4.0_dp], [3, 2])
    call assert_all_close(b, expected, 0.0_dp, 'rot90')

    v = fliplr([1.0_dp, 3.0_dp, 5.0_dp])
    call assert_all_close(v, [5.0_dp, 3.0_dp, 1.0_dp], 0.0_dp, 'vector flip')

    b = repmat(eye(2), 2, 3)
    call assert_true(all(shape(b) == [4, 6]), 'repmat shape')
    call assert_close(sum(b), 12.0_dp, 0.0_dp, 'repmat values')
    b = repmat(reshape([1.0_dp, 2.0_dp, 3.0_dp], [3, 1]), 2, 1)
    call assert_true(all(shape(b) == [2, 3]), 'source-compatible repmat shape')
    call assert_all_close(b, reshape([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp], [2, 3]), &
                          0.0_dp, 'source-compatible repmat values')
    b = reshape2d(a, 3, 2)
    call assert_true(all(shape(b) == [3, 2]), 'reshape shape')
    call assert_close(b(3, 2), 6.0_dp, 0.0_dp, 'reshape ordering')

    grid = meshgrid([1.0_dp, 2.0_dp, 3.0_dp], [10.0_dp, 20.0_dp])
    call assert_all_close(grid%x, reshape([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, &
                                           3.0_dp, 3.0_dp], [2, 3]), 0.0_dp, 'meshgrid x')
    call assert_all_close(grid%y, reshape([10.0_dp, 20.0_dp, 10.0_dp, 20.0_dp, &
                                           10.0_dp, 20.0_dp], [2, 3]), 0.0_dp, 'meshgrid y')

    b = padarray(reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2]), 1, 1, &
                 'replicate', 'both')
    call assert_true(all(shape(b) == [4, 4]), 'pad shape')
    call assert_close(b(1, 1), 1.0_dp, 0.0_dp, 'replicate pad corner')
    call assert_close(b(4, 4), 4.0_dp, 0.0_dp, 'replicate pad opposite')

    dims = shape_of(a)
    call assert_true(all(dims == [2, 3]), 'shape_of')
    call assert_int_equal(size_dim(a, 3), 1, 'singleton dimension')
    call assert_int_equal(numel(a), 6, 'numel')
    call assert_int_equal(ndims(a), 2, 'ndims')
    call assert_true(.not. isempty(a), 'isempty')
    idx = find(reshape([0.0_dp, 2.0_dp, 3.0_dp, 0.0_dp], [2, 2]))
    call assert_true(all(idx == [2, 3]), 'matrix find column-major')
    idx = find([0.0_dp, 2.0_dp, 0.0_dp, -1.0_dp])
    call assert_true(all(idx == [2, 4]), 'find')

    write(*, '(a)') 'test_array: PASS'
end program test_array
