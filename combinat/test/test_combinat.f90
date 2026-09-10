program test_combinat
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use combinat
    implicit none

    integer, allocatable :: ci(:, :)
    integer, allocatable :: c2(:, :)
    integer, allocatable :: hc_expected_i(:, :)
    integer, allocatable :: perm(:, :)
    integer, allocatable :: simplex(:, :)
    integer, allocatable :: xi(:)
    integer, allocatable :: xi_labels(:)
    integer, allocatable :: multi(:, :)
    integer, allocatable :: multz(:, :)
    real(dp), allocatable :: cr(:, :)
    real(dp), allocatable :: grid(:, :)
    real(dp), allocatable :: nv(:)
    real(dp), allocatable :: sv(:)
    real(dp), allocatable :: xr(:)
    type(combinat_rng_state) :: rng
    integer :: info

    call assert_close(fact(4.0_dp), 24.0_dp, 1.0e-12_dp, 'fact')
    call assert_close(logfact(4.0_dp), log(24.0_dp), 1.0e-12_dp, 'logfact')
    call assert_close(ncm(5.0_dp, 2.0_dp), 10.0_dp, 1.0e-12_dp, 'nCm integer')
    call assert_close(ncm(-3.0_dp, 2.0_dp), 6.0_dp, 1.0e-12_dp, 'nCm negative n')
    call assert_close(ncm(2.5_dp, 4.0_dp), -0.0390625_dp, 1.0e-12_dp, 'nCm noninteger n')
    call assert_close(ncm(3.0_dp, 4.0_dp), 0.0_dp, 1.0e-12_dp, 'nCm impossible integer choice')
    call assert_true(ieee_is_nan(ncm(5.0_dp, 2.5_dp)), 'nCm noninteger m returns NaN')

    call ncm_vec([5.0_dp, 6.0_dp], [2.0_dp], nv, info=info)
    call assert_true(info == combinat_success, 'nCm recycling status')
    call assert_real_vector(nv, [10.0_dp, 15.0_dp], 1.0e-12_dp, 'nCm recycling')

    call assert_close(dmnom([1.0_dp, 1.0_dp], [0.25_dp, 0.75_dp]), 0.375_dp, 1.0e-12_dp, 'dmnom')
    call assert_close(dmnom([1.0_dp], [0.5_dp, 0.5_dp]), 0.0_dp, 1.0e-12_dp, 'dmnom default size before recycling')

    call combn_indices(4, 2, ci, info)
    call assert_true(info == combinat_success, 'combn status')
    call assert_int_matrix(ci, reshape([1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4], [2, 6]), 'combn order')

    call combn([10.0_dp, 20.0_dp, 30.0_dp], 2, cr, info)
    call assert_true(info == combinat_success, 'combn real status')
    call assert_real_matrix(cr, reshape([10.0_dp, 20.0_dp, 10.0_dp, 30.0_dp, 20.0_dp, 30.0_dp], [2, 3]), &
        1.0e-12_dp, 'combn real values')

    call combn2_indices(4, c2, info)
    call assert_true(info == combinat_success, 'combn2 status')
    call assert_int_matrix(c2, reshape([1, 1, 1, 2, 2, 3, 2, 3, 4, 3, 4, 4], [6, 2]), 'combn2 order')

    call hcube([2, 3], grid, info=info)
    call assert_true(info == combinat_success, 'hcube status')
    allocate(hc_expected_i(6, 2))
    hc_expected_i = reshape([1, 2, 1, 2, 1, 2, 1, 1, 2, 2, 3, 3], [6, 2])
    call assert_real_matrix(grid, real(hc_expected_i, dp), 1.0e-12_dp, 'hcube lattice order')

    call hcube([2, 2], grid, scale=[10.0_dp, 100.0_dp], translation=[-1.0_dp, 0.5_dp], info=info)
    call assert_true(info == combinat_success, 'hcube affine status')
    call assert_real_matrix(grid, reshape([9.0_dp, 19.0_dp, 9.0_dp, 19.0_dp, &
        100.5_dp, 100.5_dp, 200.5_dp, 200.5_dp], [4, 2]), 1.0e-12_dp, 'hcube affine transform')

    call assert_close(nsimplex(3.0_dp, 4.0_dp), 15.0_dp, 1.0e-12_dp, 'nsimplex')
    call nsimplex_vec([3.0_dp], [2.0_dp, 3.0_dp], sv, info)
    call assert_true(info == combinat_success, 'nsimplex recycling status')
    call assert_real_vector(sv, [6.0_dp, 10.0_dp], 1.0e-12_dp, 'nsimplex recycling')

    call xsimplex(3, 2, simplex, info)
    call assert_true(info == combinat_success, 'xsimplex status')
    call assert_int_matrix(simplex, reshape([2, 0, 0, 1, 1, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 0, 2], [3, 6]), &
        'xsimplex order')
    call assert_true(all(sum(simplex, dim=1) == 2), 'xsimplex column totals')

    call permn_indices(3, perm, info)
    call assert_true(info == combinat_success, 'permn status')
    call assert_int_matrix(perm, reshape([1, 2, 3, 1, 3, 2, 3, 1, 2, 3, 2, 1, 2, 3, 1, 2, 1, 3], [3, 6]), &
        'permn minimal-change order')

    call x2u_indices([2, 0, 3], xi, info)
    call assert_true(info == combinat_success, 'x2u index status')
    call assert_int_vector(xi, [1, 1, 3, 3, 3], 'x2u indices')

    call x2u([1, 2], [10, 20], xi_labels, info)
    call assert_true(info == combinat_success, 'x2u integer-label status')
    call assert_int_vector(xi_labels, [10, 20, 20], 'x2u integer labels')

    call x2u([1, 2], [1.5_dp, 2.5_dp], xr, info)
    call assert_true(info == combinat_success, 'x2u real-label status')
    call assert_real_vector(xr, [1.5_dp, 2.5_dp, 2.5_dp], 1.0e-12_dp, 'x2u real labels')

    call rng_seed(rng, 123_int64)
    call rmultinomial([3, 2], reshape([0.2_dp, 0.5_dp, 0.3_dp, 0.25_dp, 0.5_dp, 0.25_dp], [2, 3]), &
        rng, multi, info=info)
    call assert_true(info == combinat_success, 'rmultinomial status')
    call assert_int_matrix(multi, reshape([2, 1, 0, 1, 1, 0], [2, 3]), 'rmultinomial deterministic counts')
    call assert_int_vector(sum(multi, dim=2), [3, 2], 'rmultinomial row totals')

    call rng_seed(rng, 123_int64)
    call rmultz2([3, 2], [0.2_dp, 0.3_dp, 0.5_dp], rng, multz, info=info)
    call assert_true(info == combinat_success, 'rmultz2 status')
    call assert_int_matrix(multz, reshape([2, 0, 1, 0, 1, 1], [3, 2]), 'rmultz2 deterministic counts')
    call assert_int_vector(sum(multz, dim=1), [3, 2], 'rmultz2 column totals')

    print '(a)', 'all combinat tests passed'

contains

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition !! Boolean condition that must hold for the test to pass.
        character(*), intent(in) :: label !! Human-readable test label reported on failure.

        if (.not. condition) then
            write (*, '(a)') 'FAILED: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! Computed scalar value under test.
        real(dp), intent(in) :: expected !! Reference scalar value expected from the translated R semantics.
        real(dp), intent(in) :: tolerance !! Maximum allowed absolute difference between actual and expected values.
        character(*), intent(in) :: label !! Human-readable test label reported on failure.

        call assert_true(abs(actual - expected) <= tolerance, label)
    end subroutine assert_close

    subroutine assert_int_vector(actual, expected, label)
        integer, intent(in) :: actual(:) !! Computed integer vector under test.
        integer, intent(in) :: expected(:) !! Exact reference integer vector.
        character(*), intent(in) :: label !! Human-readable test label reported on failure.

        call assert_true(size(actual) == size(expected), trim(label)//' shape')
        call assert_true(all(actual == expected), label)
    end subroutine assert_int_vector

    subroutine assert_real_vector(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual(:) !! Computed real vector under test.
        real(dp), intent(in) :: expected(:) !! Reference real vector.
        real(dp), intent(in) :: tolerance !! Maximum allowed elementwise absolute difference.
        character(*), intent(in) :: label !! Human-readable test label reported on failure.

        call assert_true(size(actual) == size(expected), trim(label)//' shape')
        call assert_true(all(abs(actual - expected) <= tolerance), label)
    end subroutine assert_real_vector

    subroutine assert_int_matrix(actual, expected, label)
        integer, intent(in) :: actual(:, :) !! Computed integer matrix under test.
        integer, intent(in) :: expected(:, :) !! Exact reference integer matrix.
        character(*), intent(in) :: label !! Human-readable test label reported on failure.

        call assert_true(all(shape(actual) == shape(expected)), trim(label)//' shape')
        call assert_true(all(actual == expected), label)
    end subroutine assert_int_matrix

    subroutine assert_real_matrix(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual(:, :) !! Computed real matrix under test.
        real(dp), intent(in) :: expected(:, :) !! Reference real matrix.
        real(dp), intent(in) :: tolerance !! Maximum allowed elementwise absolute difference.
        character(*), intent(in) :: label !! Human-readable test label reported on failure.

        call assert_true(all(shape(actual) == shape(expected)), trim(label)//' shape')
        call assert_true(all(abs(actual - expected) <= tolerance), label)
    end subroutine assert_real_matrix

end program test_combinat
